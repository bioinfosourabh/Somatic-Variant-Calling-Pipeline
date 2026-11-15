# Somatic Variant Calling Pipeline

A Bash pipeline that takes paired-end FASTQ files from tumour samples to filtered somatic variant calls (VCF) on the hg38 reference, using fastp, BWA-MEM, samtools, GATK4 (Mutect2) and bcftools.

Variants are called in **tumour-only mode**: there is no matched normal, and germline and recurrent artefactual calls are filtered using a panel of normals and a gnomAD allele-frequency germline resource. For tumour-normal calling and a containerised version, see [ngs-automation-nextflow](https://github.com/bioinfosourabh/ngs-automation-nextflow).

## Workflow

For each sample, [`scripts/somatic_variant_calling.sh`](scripts/somatic_variant_calling.sh) runs:

| Step | Tool | Output |
|---|---|---|
| Read QC and trimming | `fastp` | filtered FASTQ, HTML and JSON reports |
| Alignment and sorting | `bwa mem` with read groups, `samtools sort` | sorted, indexed BAM |
| Alignment statistics | `samtools flagstat` | `<sample>.Stat.txt` |
| Duplicate marking | `gatk MarkDuplicatesSpark` | `<sample>_markdup.bam` |
| Base quality score recalibration | `gatk BaseRecalibrator`, `gatk ApplyBQSR` (dbSNP and known-indel sites) | `<sample>_recall.bam` |
| Somatic variant calling | `gatk Mutect2` with `--germline-resource`, `--panel-of-normals` and F1R2 counts | `<sample>.vcf.gz` |
| Read-orientation artefact model | `gatk LearnReadOrientationModel` | `<sample>.read-orientation-model.tar.gz` |
| Contamination (optional, `common_variants`) | `gatk GetPileupSummaries`, `gatk CalculateContamination` | `<sample>.contamination.table` |
| Filtering | `gatk FilterMutectCalls`, then `bcftools view -f PASS` | `<sample>_filtered.vcf.gz`, `<sample>_filtered_pass.vcf.gz` |

Samples are processed one after another; each tool uses the number of threads given with `-t`.

## Requirements

- BWA, samtools, fastp, bcftools
- GATK 4.x and Java 17
- Reference files for hg38: genome FASTA with BWA and samtools indexes and a sequence dictionary, dbSNP VCF (indexed), a panel of normals VCF, and a gnomAD allele-frequency-only VCF; optionally known-indel VCFs and a common-variant VCF for contamination estimation

Example installation on Ubuntu:

```sh
sudo apt update && sudo apt install -y bwa samtools fastp bcftools openjdk-17-jre-headless
gatk_version="4.6.1.0"
wget https://github.com/broadinstitute/gatk/releases/download/$gatk_version/gatk-$gatk_version.zip
unzip gatk-$gatk_version.zip && sudo mv gatk-$gatk_version /opt/gatk
export PATH=/opt/gatk:$PATH
```

The file [`Pre-requisites`](Pre-requisites) lists download commands for the reference genome and the resource VCFs.

## Configuration

All paths are set in a configuration file, [`config.sh`](config.sh), which the script sources. Copy it, fill in the paths and pass it with `-c`:

```bash
file_ref="/path/to/reference/hg38.fa"
known_site="/path/to/vcf/dbsnp.vcf.gz"
known_indels=("/path/to/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
af_only="/path/to/vcf/af-only-gnomad.hg38.vcf.gz"
pon="/path/to/vcf/1000g_pon.hg38.vcf.gz"
common_variants=""     # e.g. small_exac_common_3.hg38.vcf.gz to estimate contamination
intervals=""           # panel/exome target BED; empty for whole genomes
```

Before anything runs, the script checks that the tools are on `PATH`, that the reference has its `.fai`, `.dict` and BWA index, and that every resource VCF exists and is indexed (`.tbi` for `.vcf.gz`, `.idx` for `.vcf`).

## Input

Paired-end, gzipped FASTQ files in one directory, named `<sample>_1.fastq.gz` and `<sample>_2.fastq.gz`. Input files are only read, never moved.

## Usage

```sh
bash scripts/somatic_variant_calling.sh -i /path/to/fastq_files -o /path/to/results -c my_config.sh -t 16
```

| Option | Description |
|---|---|
| `-i DIR` | FASTQ directory (required) |
| `-o DIR` | output directory (required) |
| `-c FILE` | configuration file (default: `config.sh` in the repository) |
| `-t N` | threads (default 8) |
| `-s LIST` | comma-separated subset of samples |
| `-n` | dry run: check inputs and list the steps that would run |
| `-f` | force: re-run all steps |

Each completed step writes a marker in `<outdir>/.steps/<sample>/`; running the same command again skips completed steps, so an interrupted run continues where it stopped. When a step re-runs, the later steps of that sample re-run too.

## Output

```
<outdir>/
├── filtered_qc_report/     fastp-filtered FASTQ files, HTML and JSON reports
├── Mapsam/                 BAM files, flagstat summaries, duplicate metrics, BQSR tables, recalibrated BAMs
├── Somatic_VCF/            raw Mutect2 VCF, orientation model, contamination tables,
│                           FilterMutectCalls VCF and PASS-only VCF per sample
├── sample_status.tsv       completed or failed, per sample and run
└── logs/                   pipeline_<timestamp>.log and one log per sample and step
```

## Notes and limitations

- Tumour-only calling is less specific than tumour-normal calling: germline variants missing from (or represented differently in) the germline resource can pass the filters. Interpret results with that in mind.
- A sample whose step fails is recorded as `failed` in `sample_status.tsv` with the end of the step log printed; the other samples continue and the script exits with status 1 at the end.
- No annotation step is included; annotate the PASS VCFs with a tool of your choice (for example VEP or ANNOVAR).
- Version 2.0 of the script moved all settings to `config.sh`, added argument parsing, `set -euo pipefail`, dependency and resource checks, per-step logs and resumable steps, the read-orientation model and optional contamination estimation, and no longer moves input FASTQ files into `Analysis_Completed/` and `Error_Samples/` (status is written to `sample_status.tsv` instead).

## License

MIT, see [LICENSE](LICENSE).

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)
