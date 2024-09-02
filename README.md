# Virus variant calling pipeline

### Usage
    git clone https://github.com/chihyuchg/virus_variant_calling.git

    conda env create -f environment.yaml
    conda activate virus_variant_calling

Run pipeline:

    nextflow run main \
    --input <input directory with fastq files>
    --output <outdir>
    --batch_name <batch_name>
    --ref_fa <reference fna> 

### Dependencies

- Python (>= 3.9)
- Gatk4 (4.3.0.0)
- hisat2 (2.2.1)
- SnpEff (5.2)

### References:
[1] https://www.ncbi.nlm.nih.gov/sra/docs/sars-cov-2-illumina-variant-calling-pipeline/
