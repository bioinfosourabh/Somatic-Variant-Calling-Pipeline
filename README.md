# Somatic Variant Calling Pipeline

## Overview
This pipeline processes paired-end FASTQ files for **somatic variant calling** using **BWA-MEM, GATK Mutect2, and bcftools**. It performs quality control, read alignment, duplicate marking, base recalibration, and variant calling, generating a **filtered VCF file** for downstream analysis.

## Features
- **End-to-end automation** for somatic variant calling
- **Quality control and filtering** for reliable variant detection
- **Parallel processing** for large-scale sequencing data
- **Logging and error handling** to track each step

## Installation
Ensure the following tools are installed:
- **BWA** (Burrows-Wheeler Aligner)
- **Samtools**
- **Fastp**
- **GATK (Genome Analysis Toolkit)**
- **bcftools**

To install dependencies (on Ubuntu):
```sh
sudo apt update && sudo apt install -y bwa samtools fastp bcftools
gatk_version="4.2.6.1"
wget https://github.com/broadinstitute/gatk/releases/download/$gatk_version/gatk-$gatk_version.zip
unzip gatk-$gatk_version.zip && mv gatk-$gatk_version /opt/gatk
export PATH=/opt/gatk:$PATH
```

## Usage
Run the pipeline by executing the script with your FASTQ directory:
```sh
bash somatic_variant_calling.sh /path/to/fastq_files
```

The script will:
1. Perform **quality control** using Fastp
2. Align reads to **hg38 reference genome** using BWA-MEM
3. Apply **duplicate marking and base recalibration** using GATK
4. Call **somatic variants** using Mutect2
5. Filter **high-confidence somatic variants**

## Output Files
- **Mapped BAM files** (with alignment statistics)
- **Filtered VCF files** (variant calls)
- **Quality control reports**
- **Log files for tracking errors and progress**

## Script
The following is the **full pipeline script**:

```bash
#!/bin/bash
# Somatic Variant Calling Pipeline
set +e  # Allow the script to continue even if some samples fail
set -o pipefail  # Ensure pipeline failures are detected

inputDir=$(pwd)
completedDir="$inputDir/Analysis_Completed"
errorDir="$inputDir/Error_Samples"
mkdir -p "$completedDir" "$errorDir"

dir_qc="$inputDir/fastqc"
dir_fastp="$inputDir/filtered_qc_report"
dir_map="$inputDir/Mapsam"
dir_vcf="$inputDir/Somatic_VCF"
mkdir -p "$dir_map" "$dir_vcf"

file_ref="/path/to/reference/hg38.fa"
known_site="/path/to/vcf/dbsnp.vcf.gz"
platform="ILLUMINA"
pon="/path/to/vcf/somatic-hg38_1000g_pon.vcf"
af_only="/path/to/vcf/somatic-hg38_af-only-gnomad.vcf"
tmp="/path/to/tmp"

timestamp=$(date +'%Y%m%d%H%M%S')
log_file="$inputDir/log_${timestamp}.txt"
terminal_log="$inputDir/terminal_log_${timestamp}.txt"
exec > >(tee -a "$terminal_log") 2>&1

echo "Script started." | tee -a "$log_file"

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
    
    fastp -i "$filep1" -I "$filep2" \
         -o "${dir_fastp}/${sample}_1_filtered.fastq.gz" \
         -O "${dir_fastp}/${sample}_2_filtered.fastq.gz" \
         --html "${dir_fastp}/${sample}.html" --thread 2
    
    bwa mem -t 32 -R "@RG\tID:${sample}\tPL:${platform}\tSM:${sample}" \
        "$file_ref" "${dir_fastp}/${sample}_1_filtered.fastq.gz" "${dir_fastp}/${sample}_2_filtered.fastq.gz" | \
        samtools sort -@ 12 -O BAM -o "${dir_map}/${sample}.bam"

    gatk AddOrReplaceReadGroups -I "${dir_map}/${sample}.bam" \
        -O "${dir_map}/${sample}_RG.bam" \
        --RGLB "Lib-${sample}" --RGPL "$platform" \
        --RGPU "${platform}_${sample}" --RGSM "$sample" \
        --VALIDATION_STRINGENCY LENIENT
  
    samtools flagstat "${dir_map}/${sample}.bam" > "${dir_map}/${sample}.Stat.txt"

    gatk MarkDuplicatesSpark -I "${dir_map}/${sample}.bam" \
        -O "${dir_map}/${sample}_markdup.bam" --spark-master local[12]
    
    gatk BaseRecalibrator -I "${dir_map}/${sample}_markdup.bam" \
        --known-sites "$known_site" -O "${dir_map}/${sample}_recall.table" -R "$file_ref"
    
    gatk ApplyBQSR --bqsr-recal-file "${dir_map}/${sample}_recall.table" \
        -I "${dir_map}/${sample}_markdup.bam" -O "${dir_map}/${sample}_recall.bam" -R "$file_ref"
    
    gatk Mutect2 -R "$file_ref" -I "${dir_map}/${sample}_recall.bam" \
        --germline-resource "$af_only" --panel-of-normals "$pon" \
        -O "${dir_vcf}/${sample}.vcf"
    
    gatk FilterMutectCalls -R "$file_ref" -V "${dir_vcf}/${sample}.vcf" \
        -O "${dir_vcf}/${sample}_filtered.vcf"
    
    bcftools view -f PASS "${dir_vcf}/${sample}_filtered.vcf" -o "${dir_vcf}/${sample}_filtered_pass.vcf"}

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
```

## License
This pipeline is released under the MIT License.

For any issues or improvements, feel free to contribute!
