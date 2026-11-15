#!/usr/bin/env bash
#
# Somatic short-variant calling, tumour-only mode (GATK Mutect2)
#
#   fastp -> BWA-MEM -> samtools sort -> MarkDuplicatesSpark -> BQSR
#   -> Mutect2 (germline resource + panel of normals, F1R2 counts)
#   -> LearnReadOrientationModel [-> GetPileupSummaries + CalculateContamination]
#   -> FilterMutectCalls -> PASS-only VCF
#
# Each finished step leaves a marker in <outdir>/.steps/<sample>/; re-running
# the command skips completed steps. A failed sample is recorded in
# <outdir>/sample_status.tsv and the other samples continue.
#
# Usage: somatic_variant_calling.sh -i FASTQ_DIR -o OUTDIR -c config.sh [-t THREADS]
#
# Author: Sourabh Kumar <bioinfosourabh@gmail.com>
set -euo pipefail

readonly VERSION="2.0.0"
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<USAGE
Somatic variant calling pipeline (tumour-only) v${VERSION}

Usage: $(basename "$0") -i FASTQ_DIR -o OUTDIR [-c CONFIG] [options]

  -i DIR     directory with paired FASTQ files <sample>_1.fastq.gz / <sample>_2.fastq.gz  (required)
  -o DIR     output directory                                                           (required)
  -c FILE    configuration file                                     [${script_dir}/../config.sh]
  -t N       threads                                                                    [8]
  -s LIST    comma-separated sample names to run (default: all samples in FASTQ_DIR)
  -n         dry run: check inputs and show the planned steps, run nothing
  -f         force: ignore step markers and re-run everything
  -h         show this help
USAGE
}

input_dir=""; outdir=""; config="${script_dir}/../config.sh"; threads=8; sample_filter=""
dry_run=false; force=false
while getopts ":i:o:c:t:s:nfh" opt; do
    case "$opt" in
        i) input_dir="$OPTARG" ;;
        o) outdir="$OPTARG" ;;
        c) config="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        s) sample_filter="$OPTARG" ;;
        n) dry_run=true ;;
        f) force=true ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: -$OPTARG needs a value" >&2; usage >&2; exit 2 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n "$input_dir" && -n "$outdir" ]] || { usage >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: -t must be a positive integer" >&2; exit 2; }
[[ -d "$input_dir" ]] || { echo "ERROR: input directory not found: $input_dir" >&2; exit 1; }
[[ -f "$config" ]] || { echo "ERROR: config file not found: $config" >&2; exit 1; }

# ------------------------------------------------------------------ config
file_ref=""; known_site=""; known_indels=(); af_only=""; pon=""; common_variants=""
platform="ILLUMINA"; intervals=""; interval_padding=100
fastp_args=(--qualified_quality_phred 20); java_mem="8g"; tmp_dir=""
# shellcheck source=/dev/null
source "$config"

input_dir=$(cd "$input_dir" && pwd)
mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)

dir_fastp="$outdir/filtered_qc_report"
dir_map="$outdir/Mapsam"
dir_vcf="$outdir/Somatic_VCF"
dir_logs="$outdir/logs"
dir_steps="$outdir/.steps"
tmp_dir="${tmp_dir:-$outdir/tmp}"
status_file="$outdir/sample_status.tsv"
mkdir -p "$dir_fastp" "$dir_map" "$dir_vcf" "$dir_logs" "$dir_steps" "$tmp_dir"

main_log="$dir_logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" | tee -a "$main_log" >&2; }
die() { log "ERROR: $*"; exit 1; }

# ------------------------------------------------------------------ checks
check_dependencies() {
    local missing=() tool
    for tool in fastp bwa samtools gatk bcftools java; do
        command -v "$tool" >/dev/null 2>&1 || missing+=("$tool")
    done
    if (( ${#missing[@]} > 0 )); then
        if [[ "$dry_run" == true ]]; then
            log "WARNING: not on PATH: ${missing[*]}"
        else
            die "required tools not on PATH: ${missing[*]}"
        fi
    fi
}

vcf_index() {
    if [[ "$1" == *.gz ]]; then echo "$1.tbi"; else echo "$1.idx"; fi
}

require_vcf() {
    local label="$1" path="$2"
    [[ -n "$path" ]] || die "$label is not set in $config"
    [[ -f "$path" ]] || die "$label not found: $path"
    [[ -f "$(vcf_index "$path")" ]] || die "$label has no index ($(vcf_index "$path"))"
}

check_resources() {
    [[ -n "$file_ref" && -f "$file_ref" ]] || die "file_ref not found: '$file_ref'"
    [[ -f "${file_ref}.fai" ]] || die "missing ${file_ref}.fai (samtools faidx)"
    [[ -f "${file_ref%.*}.dict" ]] || die "missing ${file_ref%.*}.dict (gatk CreateSequenceDictionary)"
    [[ -f "${file_ref}.bwt" ]] || die "missing BWA index for $file_ref (bwa index)"
    require_vcf known_site "$known_site"
    local v
    for v in "${known_indels[@]}"; do require_vcf known_indels "$v"; done
    require_vcf af_only "$af_only"
    if [[ -n "$pon" ]]; then
        require_vcf pon "$pon"
    else
        log "WARNING: no panel of normals (pon); expect more recurrent artefacts"
    fi
    [[ -z "$common_variants" ]] || require_vcf common_variants "$common_variants"
    [[ -z "$intervals" || -f "$intervals" ]] || die "intervals file not found: $intervals"
}

# ------------------------------------------------------------------ step runner
declare -A scope_rerun=()
run_step() {
    local name="$1" scope="$2"; shift 2
    local marker="$dir_steps/$scope/$name.done" step_log="$dir_logs/$scope/$name.log"
    mkdir -p "$dir_steps/$scope" "$dir_logs/$scope"
    if [[ -f "$marker" && "$force" == false && -z "${scope_rerun[$scope]:-}" ]]; then
        log "$scope: skip $name (done)"
        return 0
    fi
    if [[ "$dry_run" == true ]]; then
        log "$scope: would run $name"
        return 0
    fi
    scope_rerun[$scope]=1
    rm -f "$marker"
    log "$scope: $name"
    if "$@" > "$step_log" 2>&1; then
        touch "$marker"
    else
        log "$scope: $name FAILED; last lines of $step_log:"
        tail -n 8 "$step_log" | sed 's/^/    /' | tee -a "$main_log" >&2
        return 1
    fi
}

gatk_run() {
    gatk --java-options "-Xmx${java_mem} -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${tmp_dir}" "$@"
}

interval_args=()
if [[ -n "$intervals" ]]; then
    interval_args=(-L "$intervals" --interval-padding "$interval_padding")
fi

# ------------------------------------------------------------------ steps
# Commands are chained with && so that any failure is returned to run_step.
step_fastp() {
    local s="$1"
    fastp -i "$input_dir/${s}_1.fastq.gz" -I "$input_dir/${s}_2.fastq.gz" \
        -o "$dir_fastp/${s}_1_filtered.fastq.gz" -O "$dir_fastp/${s}_2_filtered.fastq.gz" \
        --html "$dir_fastp/${s}.html" --json "$dir_fastp/${s}.json" \
        --thread "$(( threads < 16 ? threads : 16 ))" "${fastp_args[@]}"
}

step_align() {
    local s="$1"
    local bam="$dir_map/${s}.bam"
    bwa mem -t "$threads" \
        -R "@RG\tID:${s}\tPL:${platform}\tLB:LIB_${s}\tPU:${platform}_${s}\tSM:${s}" \
        "$file_ref" "$dir_fastp/${s}_1_filtered.fastq.gz" "$dir_fastp/${s}_2_filtered.fastq.gz" \
        | samtools sort -@ "$threads" -T "$tmp_dir/${s}.sort" -O BAM -o "$bam" - &&
    samtools index "$bam" &&
    samtools flagstat "$bam" > "$dir_map/${s}.Stat.txt" &&
    if (( $(samtools view -c -F 0x904 "$bam") == 0 )); then
        echo "no mapped reads in $bam" >&2
        return 1
    fi
}

step_markdup() {
    local s="$1"
    gatk_run MarkDuplicatesSpark -I "$dir_map/${s}.bam" -O "$dir_map/${s}_markdup.bam" \
        -M "$dir_map/${s}_markdup_metrics.txt" \
        --spark-master "local[${threads}]" --tmp-dir "$tmp_dir"
}

step_bqsr() {
    local s="$1" ks=(--known-sites "$known_site") v
    for v in "${known_indels[@]}"; do ks+=(--known-sites "$v"); done
    gatk_run BaseRecalibrator -R "$file_ref" -I "$dir_map/${s}_markdup.bam" "${ks[@]}" \
        "${interval_args[@]}" -O "$dir_map/${s}_recall.table" &&
    gatk_run ApplyBQSR -R "$file_ref" -I "$dir_map/${s}_markdup.bam" \
        --bqsr-recal-file "$dir_map/${s}_recall.table" -O "$dir_map/${s}_recall.bam"
}

step_mutect2() {
    local s="$1" pon_args=()
    [[ -z "$pon" ]] || pon_args=(--panel-of-normals "$pon")
    gatk_run Mutect2 -R "$file_ref" -I "$dir_map/${s}_recall.bam" \
        --germline-resource "$af_only" "${pon_args[@]}" "${interval_args[@]}" \
        --f1r2-tar-gz "$dir_vcf/${s}.f1r2.tar.gz" \
        --native-pair-hmm-threads "$threads" \
        -O "$dir_vcf/${s}.vcf.gz" &&
    gatk_run LearnReadOrientationModel -I "$dir_vcf/${s}.f1r2.tar.gz" \
        -O "$dir_vcf/${s}.read-orientation-model.tar.gz"
}

step_contamination() {
    local s="$1" iv=(-L "$common_variants")
    [[ -z "$intervals" ]] || iv+=(-L "$intervals" --interval-set-rule INTERSECTION)
    gatk_run GetPileupSummaries -I "$dir_map/${s}_recall.bam" -V "$common_variants" "${iv[@]}" \
        -O "$dir_vcf/${s}.pileups.table" &&
    gatk_run CalculateContamination -I "$dir_vcf/${s}.pileups.table" \
        --tumor-segmentation "$dir_vcf/${s}.segments.table" \
        -O "$dir_vcf/${s}.contamination.table"
}

step_filter() {
    local s="$1" extra=()
    if [[ -n "$common_variants" ]]; then
        extra=(--contamination-table "$dir_vcf/${s}.contamination.table"
               --tumor-segmentation "$dir_vcf/${s}.segments.table")
    fi
    gatk_run FilterMutectCalls -R "$file_ref" -V "$dir_vcf/${s}.vcf.gz" \
        --ob-priors "$dir_vcf/${s}.read-orientation-model.tar.gz" "${extra[@]}" \
        -O "$dir_vcf/${s}_filtered.vcf.gz" &&
    bcftools view -f PASS -Oz -o "$dir_vcf/${s}_filtered_pass.vcf.gz" "$dir_vcf/${s}_filtered.vcf.gz" &&
    bcftools index -f -t "$dir_vcf/${s}_filtered_pass.vcf.gz"
}

process_sample() {
    local s="$1"
    run_step fastp "$s" step_fastp "$s" &&
    run_step align "$s" step_align "$s" &&
    run_step markdup "$s" step_markdup "$s" &&
    run_step bqsr "$s" step_bqsr "$s" &&
    run_step mutect2 "$s" step_mutect2 "$s" &&
    if [[ -n "$common_variants" ]]; then run_step contamination "$s" step_contamination "$s"; fi &&
    run_step filter "$s" step_filter "$s"
}

# ------------------------------------------------------------------ main
log "somatic-variant-calling v${VERSION} (tumour-only)"
log "input: $input_dir  output: $outdir  config: $config  threads: $threads"
check_dependencies
check_resources

samples=()
for f in "$input_dir"/*_1.fastq.gz; do
    [[ -e "$f" ]] || continue
    s=$(basename "$f" _1.fastq.gz)
    if [[ ! -f "$input_dir/${s}_2.fastq.gz" ]]; then
        log "WARNING: skipping $s: ${s}_2.fastq.gz not found"
        continue
    fi
    if [[ -n "$sample_filter" && ",${sample_filter}," != *",${s},"* ]]; then
        continue
    fi
    samples+=("$s")
done
(( ${#samples[@]} > 0 )) || die "no FASTQ pairs (<sample>_1.fastq.gz, <sample>_2.fastq.gz) found in $input_dir"
log "samples (${#samples[@]}): ${samples[*]}"

failed=()
[[ -f "$status_file" ]] || printf 'sample\tstatus\tfinished\n' > "$status_file"
for s in "${samples[@]}"; do
    status=0
    process_sample "$s" || status=$?
    [[ "$dry_run" == true ]] && continue
    if (( status == 0 )); then
        printf '%s\tcompleted\t%s\n' "$s" "$(date '+%F %T')" >> "$status_file"
    else
        failed+=("$s")
        printf '%s\tfailed\t%s\n' "$s" "$(date '+%F %T')" >> "$status_file"
        log "$s: failed; continuing with the other samples"
    fi
done

if [[ "$dry_run" == true ]]; then
    log "dry run complete: nothing was run"
    exit 0
fi
if (( ${#failed[@]} > 0 )); then
    log "finished with failed sample(s): ${failed[*]} (see $dir_logs/<sample>/)"
    exit 1
fi
log "finished: ${#samples[@]} sample(s); PASS calls in $dir_vcf/<sample>_filtered_pass.vcf.gz"
