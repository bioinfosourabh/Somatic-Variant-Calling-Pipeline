# Somatic Variant Calling Pipeline

A Bash pipeline that takes paired-end FASTQ files from tumour samples to filtered somatic variant calls (VCF) on the hg38 reference, using fastp, BWA-MEM, samtools, GATK4 (Mutect2) and bcftools.

Variants are called in **tumour-only mode**: there is no matched normal, and germline and recurrent artefactual calls are filtered using a panel of normals and a gnomAD allele-frequency germline resource.

## Workflow

For each sample, [`scripts/somatic_variant_calling.sh`](scripts/somatic_variant_calling.sh) runs:

| Step | Tool | Output |
|---|---|---|
| Read QC and trimming | `fastp` | filtered FASTQ and HTML report |
| Alignment and sorting | `bwa mem` with read groups, `samtools sort` | sorted BAM |
| Alignment statistics | `samtools flagstat` | `<sample>.Stat.txt` |
| Duplicate marking | `gatk MarkDuplicatesSpark` | `<sample>_markdup.bam` |
| Base quality score recalibration | `gatk BaseRecalibrator`, `gatk ApplyBQSR` (dbSNP known sites) | `<sample>_recall.bam` |
| Somatic variant calling | `gatk Mutect2` with `--panel-of-normals` and `--germline-resource` | `<sample>.vcf` |
| Filtering | `gatk FilterMutectCalls`, then `bcftools view -f PASS` | `<sample>_filtered.vcf`, `<sample>_filtered_pass.vcf` |

Samples are processed one after another; individual tools use multiple threads (BWA 32, samtools 12, Spark 12 local cores). Edit these values in the script to match your machine.

## Requirements

- BWA, samtools, fastp, bcftools
- GATK 4.x
- Reference files for hg38: genome FASTA with BWA and samtools indexes and a sequence dictionary, dbSNP VCF (indexed), a panel of normals VCF, and a gnomAD allele-frequency-only VCF

Example installation on Ubuntu:

```sh
sudo apt update && sudo apt install -y bwa samtools fastp bcftools
gatk_version="4.2.6.1"
wget https://github.com/broadinstitute/gatk/releases/download/$gatk_version/gatk-$gatk_version.zip
unzip gatk-$gatk_version.zip && mv gatk-$gatk_version /opt/gatk
export PATH=/opt/gatk:$PATH
```

The file [`Pre-requisites`](Pre-requisites) lists download commands for the reference genome and the resource VCFs.

## Configuration

Set the paths at the top of the script before running:

```bash
file_ref="/path/to/reference/hg38.fa"
known_site="/path/to/vcf/dbsnp.vcf.gz"
pon="/path/to/vcf/somatic-hg38_1000g_pon.vcf"
af_only="/path/to/vcf/somatic-hg38_af-only-gnomad.vcf"
```

## Input

Paired-end, gzipped FASTQ files in one directory, named `<sample>_1.fastq.gz` and `<sample>_2.fastq.gz`.

## Usage

The script works on the **current working directory**, so run it from the folder that contains the FASTQ files:

```sh
cd /path/to/fastq_files
bash /path/to/somatic_variant_calling.sh
```

## Output

Created inside the FASTQ directory:

```
filtered_qc_report/     fastp-filtered FASTQ files and HTML QC reports
Mapsam/                 BAM files, flagstat summaries, BQSR tables, recalibrated BAMs
Somatic_VCF/            raw Mutect2 VCF, FilterMutectCalls VCF, PASS-only VCF per sample
Analysis_Completed/     input FASTQ pairs of samples that finished
Error_Samples/          input FASTQ pairs of samples that failed
log_<timestamp>.txt, terminal_log_<timestamp>.txt
```

## Notes and limitations

- Tumour-only calling is less specific than tumour-normal calling; results should be interpreted with that in mind.
- The script does not stop on errors (`set +e`). A sample is moved to `Error_Samples/` if its FASTQ pair is missing or the final filtering step fails; check the logs for errors in intermediate steps.
- No annotation step is included; annotate the PASS VCFs with a tool of your choice (for example VEP or ANNOVAR).

## License

MIT, see [LICENSE](LICENSE).

## Author

Sourabh Kumar · [bioinfosourabh@gmail.com](mailto:bioinfosourabh@gmail.com) · [bioinfosourabh.netlify.app](https://bioinfosourabh.netlify.app) · ORCID [0009-0009-5840-9830](https://orcid.org/0009-0009-5840-9830)
