#!/bin/bash
# Somatic Variant Calling Pipeline
# This script processes paired-end FASTQ files for somatic variant calling using BWA, GATK, and bcftools.

set +e  # Allow the script to continue even if some samples fail
set -o pipefail  # Ensure pipeline failures are detected

#############################
# Define Directories & Files
#############################

inputDir=$(pwd)
completedDir="$inputDir/Analysis_Completed"
errorDir="$inputDir/Error_Samples"
mkdir -p "$completedDir" "$errorDir"

dir_qc="$inputDir/fastqc"
dir_fastp="$inputDir/filtered_qc_report"
dir_map="$inputDir/Mapsam"
dir_vcf="$inputDir/Somatic_VCF"
mkdir -p "$dir_map" "$dir_vcf"

# Reference Files
file_ref="/path/to/reference/hg38.fa"
known_site="/path/to/vcf/dbsnp.vcf.gz"
platform="ILLUMINA"
pon="/path/to/vcf/somatic-hg38_1000g_pon.vcf"
af_only="/path/to/vcf/somatic-hg38_af-only-gnomad.vcf"
tmp="/path/to/tmp"

# Logging Setup
timestamp=$(date +'%Y%m%d%H%M%S')
log_file="$inputDir/log_${timestamp}.txt"
terminal_log="$inputDir/terminal_log_${timestamp}.txt"
exec > >(tee -a "$terminal_log") 2>&1

echo "Script started." | tee -a "$log_file"

##########################################
# Discover Unique Samples from FASTQ Files
##########################################

samples=()
for f in "$inputDir"/*_1.fastq.gz; do
    if [ -f "$f" ]; then
        sampleName=$(basename "$f" "_1.fastq.gz")
        samples+=("$sampleName")
    fi
done

unique_samples=($(printf "%s\n" "${samples[@]}" | sort -u))
numSamples=${#unique_samples[@]}
echo "Total samples: $numSamples" | tee -a "$log_file"

############################################
# Function: process_sample
############################################
process_sample() {
    local sample="$1"
    echo "Processing sample: $sample" | tee -a "$log_file"
    
    local filep1="$inputDir/${sample}_1.fastq.gz"
    local filep2="$inputDir/${sample}_2.fastq.gz"
    
    if [ ! -f "$filep1" ] || [ ! -f "$filep2" ]; then
        echo "ERROR: Missing files for $sample" | tee -a "$log_file"
        return 1
    fi
    
    mkdir -p "$dir_qc" "$dir_fastp"
    
    ##################
    # Quality Control
    ##################
    fastp -i "$filep1" -I "$filep2" \
         -o "${dir_fastp}/${sample}_1_filtered.fastq.gz" \
         -O "${dir_fastp}/${sample}_2_filtered.fastq.gz" \
         --html "${dir_fastp}/${sample}.html" --thread 2
    
    ##################
    # Alignment & Sorting
    ##################
    bwa mem -t 32 -R "@RG\tID:${sample}\tPL:${platform}\tSM:${sample}" \
        "$file_ref" "${dir_fastp}/${sample}_1_filtered.fastq.gz" "${dir_fastp}/${sample}_2_filtered.fastq.gz" | \
        samtools sort -@ 12 -O BAM -o "${dir_map}/${sample}.bam"
    
    samtools flagstat "${dir_map}/${sample}.bam" > "${dir_map}/${sample}.Stat.txt"
    
    ###############################
    # Post-alignment Processing
    ###############################
    gatk MarkDuplicatesSpark -I "${dir_map}/${sample}.bam" \
        -O "${dir_map}/${sample}_markdup.bam" --spark-master local[12]
    
    gatk BaseRecalibrator -I "${dir_map}/${sample}_markdup.bam" \
        --known-sites "$known_site" -O "${dir_map}/${sample}_recall.table" -R "$file_ref"
    
    gatk ApplyBQSR --bqsr-recal-file "${dir_map}/${sample}_recall.table" \
        -I "${dir_map}/${sample}_markdup.bam" -O "${dir_map}/${sample}_recall.bam" -R "$file_ref"
    
    ##################
    # Variant Calling
    ##################
    gatk Mutect2 -R "$file_ref" -I "${dir_map}/${sample}_recall.bam" \
        --germline-resource "$af_only" --panel-of-normals "$pon" \
        -O "${dir_vcf}/${sample}.vcf"
    
    gatk FilterMutectCalls -R "$file_ref" -V "${dir_vcf}/${sample}.vcf" \
        -O "${dir_vcf}/${sample}_filtered.vcf"
    
    bcftools view -f PASS "${dir_vcf}/${sample}_filtered.vcf" -o "${dir_vcf}/${sample}_filtered_pass.vcf"
}

######################################
# Main Loop: Process Each Sample
######################################
for sample in "${unique_samples[@]}"; do
    if process_sample "$sample"; then
        mv "$inputDir/${sample}_1.fastq.gz" "$completedDir/"
        mv "$inputDir/${sample}_2.fastq.gz" "$completedDir/"
    else
        mv "$inputDir/${sample}_1.fastq.gz" "$errorDir/"
        mv "$inputDir/${sample}_2.fastq.gz" "$errorDir/"
    fi
done

echo "Pipeline completed." | tee -a "$log_file"
