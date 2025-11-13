import pandas as pd
import os

configfile: "config/config.yaml"

# Load metadata
METADATA = pd.read_csv(config["metadata_file"].strip())
SAMPLES = METADATA["sample_id"].tolist()

# Map sample IDs to raw reads
READS = {row.sample_id: {"r1": row.fwd_reads, "r2": row.rev_reads} for idx, row in METADATA.iterrows()}

# Rule all: final outputs
rule all:
    input:
        # Host-filtered reads
        expand(os.path.join(config["output_dir"], "nonhost/{sample}_R1_clean.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(config["output_dir"], "nonhost/{sample}_R2_clean.fastq.gz"), sample=SAMPLES),
        # FastQC on non-host reads
        expand(os.path.join(config["output_dir"], "FastQC/{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(os.path.join(config["output_dir"], "FastQC/{sample}_R2_fastqc.html"), sample=SAMPLES),
        # MultiQC report
        [os.path.join(config["output_dir"], "FASTQC_reports", "multiqc_report.html")]

# Trimming reads
rule trim:
    input:
        r1=lambda wildcards: READS[wildcards.sample]["r1"],
        r2=lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        r1_trimmed=os.path.join(config["output_dir"], "trimmed/{sample}_R1_trimmed.fastq.gz"),
        r2_trimmed=os.path.join(config["output_dir"], "trimmed/{sample}_R2_trimmed.fastq.gz"),
        summary=os.path.join(config["output_dir"], "trimmed/{sample}_trimmomatic_summary.txt")
    conda:
        "envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
        {output.r1_trimmed} /dev/null \
        {output.r2_trimmed} /dev/null \
        -summary {output.summary} \
        ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:50 SLIDINGWINDOW:4:20
        """


#FastQC is used for quality control
rule fastQC:
    input:
        r1=os.path.join(config["output_dir"], "trimmed/{sample}_R1_trimmed.fastq.gz"),
        r2=os.path.join(config["output_dir"], "trimmed/{sample}_R2_trimmed.fastq.gz")
    output:
        fastqc_output_dir=directory(os.path.join(config["output_dir"], "FastQC_reports/{sample}")),
        r1_html=os.path.join(config["output_dir"], "FastQC/{sample}_R1_fastqc.html"),
        r2_html=os.path.join(config["output_dir"], "FastQC/{sample}_R2_fastqc.html")
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {output.fastqc_output_dir}
        fastqc {input.r1} {input.r2} -o {output.fastqc_output_dir}
        """

#MultiqC compiles results
rule multiQC:
    input:
        fastqc_output_dir=expand(os.path.join(config["output_dir"], "FastQC_reports/{sample}"),sample=SAMPLES) 
    output:
        multiqc_report_name=[os.path.join(config["output_dir"], "FASTQC_reports", "multiqc_report.html")]
    params:
        fastqc_search_dir=directory(os.path.join(config["output_dir"], "FastQC_reports"))
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        multiqc {params.fastqc_search_dir} --outdir {input.fastqc_output_dir} --filename {output.multiqc_report_name}  --dirs-depth 2
        """

# Host filtering with Hostile
rule hostile:
    input:
        r1=os.path.join(config["output_dir"], "trimmed/{sample}_R1_trimmed.fastq.gz"),
        r2=os.path.join(config["output_dir"], "trimmed/{sample}_R2_trimmed.fastq.gz")
    output:
        r1_nonhost=os.path.join(config["output_dir"], "nonhost/{sample}_R1_clean.fastq.gz"),
        r2_nonhost=os.path.join(config["output_dir"], "nonhost/{sample}_R2_clean.fastq.gz"),
        report=os.path.join(config["output_dir"], "nonhost/{sample}_hostile_report.html")
    params:
        hostile_out_dir=directory(os.path.join(config["output_dir"], "nonhost")),
        index=config["hostile_dir"],
        threads=8
    conda:
        "envs/hostile.yaml"
    shell:
        """
        hostile clean \
            --fastq1 {input.r1} --fastq2 {input.r2} \
            -o {params.hostile_out_dir} \
            --threads {params.threads} \
            --index {params.index} \
            > {output.report}
        """

