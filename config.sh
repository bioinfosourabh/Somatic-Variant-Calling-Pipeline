#!/usr/bin/env bash
# Configuration for scripts/somatic_variant_calling.sh
# Copy this file, fill in the paths for your system, and pass it with -c.
# shellcheck disable=SC2034  # variables are used by the sourcing script

# ---------------------------------------------------------------- reference
# FASTA indexed with: bwa index, samtools faidx, gatk CreateSequenceDictionary
file_ref="/path/to/reference/hg38.fa"

# ---------------------------------------------------------------- BQSR known sites
# VCFs need an index: .tbi for .vcf.gz, .idx for plain .vcf
known_site="/path/to/resources/Homo_sapiens_assembly38.dbsnp138.vcf"
known_indels=(
    "/path/to/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
)

# ---------------------------------------------------------------- Mutect2 resources
af_only="/path/to/resources/af-only-gnomad.hg38.vcf.gz"      # germline resource (required)
pon="/path/to/resources/1000g_pon.hg38.vcf.gz"               # panel of normals (recommended; empty to skip)
common_variants=""                                            # e.g. small_exac_common_3.hg38.vcf.gz; enables contamination estimation

# ---------------------------------------------------------------- data
platform="ILLUMINA"
intervals=""            # target BED for panels/exomes; empty for whole genomes
interval_padding=100

# ---------------------------------------------------------------- tools
fastp_args=(--qualified_quality_phred 20)
java_mem="8g"
tmp_dir=""              # scratch for GATK/Spark; default: <outdir>/tmp
