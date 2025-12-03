# Snakefile: Shotgun Metagenomics Pipeline with Pre-/Post-Trim FastQC, MultiQC, 
# Host Filtering, Trimming, Kraken2, Bracken, and Aggregation 
import os 
import pandas as pd 

configfile: "config/config.yaml" 

# Load metadata 
METADATA = pd.read_csv(config["metadata_file"].strip()) 
SAMPLES = METADATA["sample_id"].tolist() 

# Map sample IDs to raw reads 
READS = { 
    row.sample_id: { 
        "r1": row.fwd_reads, 
        "r2": row.rev_reads 
    } 
    for idx, row in METADATA.iterrows() 
} 

# Convenience 
pj = os.path.join 

# Directories 
trim_trunc_path = pj(config["output_dir"], "trimmed") 
nonhost_path = pj(config["output_dir"], "nonhost") 
kraken_path = pj(config["output_dir"], "kraken") 
fastqc_path = pj(config["output_dir"], "FastQC_reports") 
kraken_db_loc = config["kraken_db_loc"] 

# ------------------------------------------------------ 
# RULE ALL 
# ------------------------------------------------------ 
rule all: 
    input: 
        # Pre-trim FastQC 
        expand(pj(fastqc_path, "pretrim/{sample}", "{sample}_R1_fastqc.html"), sample=SAMPLES), 
        expand(pj(fastqc_path, "pretrim/{sample}", "{sample}_R2_fastqc.html"), sample=SAMPLES), 
        
        # MultiQC pre-trim 
        pj(fastqc_path, "pretrim_multiqc_report.html"), 
        
        # Unzipped reads 
        expand(pj(config["output_dir"], "unzipped/{sample}_R1.fastq"), sample=SAMPLES), 
        expand(pj(config["output_dir"], "unzipped/{sample}_R2.fastq"), sample=SAMPLES), 
        
        # Trimmed reads 
        expand(pj(trim_trunc_path, "{sample}_R1_trimmed.fastq"), sample=SAMPLES), 
        expand(pj(trim_trunc_path, "{sample}_R2_trimmed.fastq"), sample=SAMPLES), 
        
        # Host-filtered 
        expand(pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"), sample=SAMPLES), 
        expand(pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz"), sample=SAMPLES), 
        
        # Post-trim FastQC 
        expand(pj(fastqc_path, "{sample}", "{sample}_R1_trimmed_fastqc.html"), sample=SAMPLES), 
        expand(pj(fastqc_path, "{sample}", "{sample}_R2_trimmed_fastqc.html"), sample=SAMPLES), 
        
        # MultiQC post-trim pj(fastqc_path, "multiqc_report.html"), 
        
        # Kraken2 + Bracken 
        expand(pj(kraken_path, "{sample}.kraken.txt"), sample=SAMPLES), 
        expand(pj(kraken_path, "{sample}.kreport2"), sample=SAMPLES), 
        expand(pj(kraken_path, "{sample}.bracken.tsv"), sample=SAMPLES), 
        
        # Aggregated taxonomy output pj(kraken_path, "Combined-taxonomy.tsv") 
        
# ------------------------------------------------------ 
# UNZIP 
# ------------------------------------------------------ 
rule unzip: 
    input:
        r1=lambda wc: READS[wc.sample]["r1"], 
        r2=lambda wc: READS[wc.sample]["r2"]
    output:
        r1=pj(config["output_dir"], "unzipped/{sample}_R1.fastq"), 
        r2=pj(config["output_dir"], "unzipped/{sample}_R2.fastq") 
    shell:
        """ 
        mkdir -p {config[output_dir]}/unzipped 
        gunzip -c {input.r1} > {output.r1} 
        gunzip -c {input.r2} > {output.r2} 
        """ 
                
# ------------------------------------------------------ 
# PRE-TRIM FASTQC 
# ------------------------------------------------------ 
rule fastqc_pretrim:
    input:
        r1=pj(config["output_dir"], "unzipped/{sample}_R1.fastq"), 
        r2=pj(config["output_dir"], "unzipped/{sample}_R2.fastq") 
    output:
        dir=directory(pj(fastqc_path, "pretrim/{sample}")), 
        r1_html=pj(fastqc_path, "pretrim/{sample}", "{sample}_R1_fastqc.html"), 
        r2_html=pj(fastqc_path, "pretrim/{sample}", "{sample}_R2_fastqc.html") 
    conda: "envs/fastqc.yaml" 
    shell:
        """ 
        mkdir -p {output.dir} 
        fastqc {input.r1} {input.r2} -o {output.dir} 
        """ 
# ------------------------------------------------------ 
# MULTIQC PRETRIM 
# ------------------------------------------------------ 
rule multiqc_pretrim: 
    input: 
        expand(pj(fastqc_path, "pretrim/{sample}"), sample=SAMPLES) 
    output: 
        html=pj(fastqc_path, "pretrim_multiqc_report.html") 
    conda: "envs/fastqc.yaml" 
    shell:
        """ 
        mkdir -p $(dirname {output.html}) 
        multiqc {fastqc_path}/pretrim -o $(dirname {output.html}) 
        mv $(dirname {output.html})/multiqc_report.html {output.html} 
        """ 
# ------------------------------------------------------ 
# TRIMMING (BBDUK) 
# ------------------------------------------------------ 
rule run_bbduk:
    input: 
        r1=pj(config["output_dir"], "unzipped/{sample}_R1.fastq"), 
        r2=pj(config["output_dir"], "unzipped/{sample}_R2.fastq") 
    output: 
        r1_trimmed=pj(trim_trunc_path, "{sample}_R1_trimmed.fastq"), 
        r2_trimmed=pj(trim_trunc_path, "{sample}_R2_trimmed.fastq") 
    conda: "envs/bbmap.yaml" 
    params: 
        adapter_refs=config["adapter_dir"], 
        trim_left=config["trim_left"],
        trim_right=config["trim_right"], 
        kmer_length=config["kmer_length"] 
    shell: 
        """ 
        mkdir -p {trim_trunc_path} 
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={params.adapter_refs} \
            out1={output.r1_trimmed} out2={output.r2_trimmed} \
            ktrim=lr k={params.kmer_length} mink=4 minlen=20 qtrim=f \
            ftl={params.trim_left} ftr={params.trim_right} threads=14 overwrite=true 
        """ 
# ------------------------------------------------------ 
# HOST FILTERING (HOSTILE) 
# ------------------------------------------------------ 
rule hostile: 
    input: 
        r1=pj(trim_trunc_path, "{sample}_R1_trimmed.fastq"), 
        r2=pj(trim_trunc_path, "{sample}_R2_trimmed.fastq") 
    output: 
        r1_nonhost=pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"), 
        r2_nonhost=pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz"), 
        report=pj(nonhost_path, "{sample}_hostile_report.html") 
    params:
        outdir=directory(nonhost_path),
        index=config["hostile_dir"],
        threads=8 
    conda: "envs/hostile.yaml" 
    shell:
        """
        hostile clean --fastq1 {input.r1} --fastq2 {input.r2} \
            -o {params.outdir} --threads {params.threads} --index {params.index} \
            > {output.report} 
        """ 

# ------------------------------------------------------ 
# POST-TRIM FASTQC 
# ------------------------------------------------------ 
rule fastqc_posttrim: 
    input: 
        r1=pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"), 
        r2=pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz") 
    output: 
        dir=directory(pj(fastqc_path, "{sample}")), 
        r1_html=pj(fastqc_path, "{sample}", "{sample}_R1_trimmed_fastqc.html"), 
        r2_html=pj(fastqc_path, "{sample}", "{sample}_R2_trimmed_fastqc.html") 
    conda: "envs/fastqc.yaml" 
    shell: 
        """
        mkdir -p {output.dir}
        fastqc {input.r1} {input.r2} -o {output.dir}
        """ 

# ------------------------------------------------------ 
# MULTIQC POST-TRIM 
# ------------------------------------------------------ 
rule multiqc_posttrim: 
    input: 
        expand(pj(fastqc_path, "{sample}"), sample=SAMPLES) 
    output: 
        html=pj(fastqc_path, "multiqc_report.html") 
    conda: "envs/fastqc.yaml"
    shell: 
        """
        mkdir -p $(dirname {output.html}) 
        multiqc {fastqc_path} -o $(dirname {output.html}) 
        mv $(dirname {output.html})/multiqc_report.html {output.html} 
        """ 

# ------------------------------------------------------ 
# KRAKEN2 for taxonomic classification 
# ------------------------------------------------------ 
rule run_kraken: 
    input:
        FWD=pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"), 
        REV=pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz"), 
        HASH=pj(kraken_db_loc, "hash.k2d") 
    output:
        OUTFILE=pj(kraken_path, "{sample}.kraken.txt"),
        REPORT=pj(kraken_path, "{sample}.kreport2") 
    conda: "envs/kraken.yaml" 
    params:
        out_dir=kraken_path,
        database=kraken_db_loc 
    threads: 16
    shell: 
        """ 
        mkdir -p {params.out_dir} 
        kraken2 \
            --gzip-compressed \
            --paired \
            --db {params.database} \
            --threads {threads} \
            --output {output.OUTFILE} \
            --report {output.REPORT} \
            {input.FWD} {input.REV} 
        """ 
# ------------------------------------------------------ 
# Bracken for abundance estimation 
# ------------------------------------------------------ 
rule bracken: 
    input: 
        kraken_report=pj(kraken_path, "{sample}.kreport2")
    output: 
        bracken_report=pj(kraken_path, "{sample}.bracken.tsv") 
    params:
        level="S",
        reads=100000 
    threads: 4 
    conda: "envs/kraken.yaml" 
    shell: 
        """ 
        bracken \
            -d {config[kraken_db_loc]} \
            -i {input.kraken_report} \
            -o {output.bracken_report} \
            -l {params.level} \
            -r {params.reads} 
        """ 
# ------------------------------------------------------ 
# Bracken Aggregation 
# ------------------------------------------------------ 
rule aggregate_bracken: 
    input: 
        expand(pj(kraken_path, "{sample}.bracken.tsv"), sample=SAMPLES) 
    output: 
        html=pj(kraken_path, "Combined-taxonomy.tsv")
    conda: "envs/kraken.yaml" 
    shell: 
        """
        echo -e "sample name taxid abs frac" > {output}
        for f in {input}; do \
            s=$(basename $f .bracken.tsv); \
            awk -v s=$s 'NR>1{print s" "$0}' $f >> {output}; \
        done 
        """