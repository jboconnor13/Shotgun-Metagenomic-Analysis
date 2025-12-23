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
humann_path = pj(config["output_dir"], "humann")

# ------------------------------------------------------
# RULE ALL
# ------------------------------------------------------
rule all:
    input:
        # Unzipped reads
        expand(pj(config["output_dir"], "unzipped/{sample}_R1.fastq"), sample=SAMPLES),
        expand(pj(config["output_dir"], "unzipped/{sample}_R2.fastq"), sample=SAMPLES),
        
        # Pre-trim FastQC
        expand(pj(fastqc_path, "pretrim/{sample}", "{sample}_R1_fastqc.html"), sample=SAMPLES),
        expand(pj(fastqc_path, "pretrim/{sample}", "{sample}_R2_fastqc.html"), sample=SAMPLES),

        # MultiQC pre-trim
        pj(fastqc_path, "pretrim_multiqc_report.html"),

        # Trimmed reads
        expand(pj(trim_trunc_path, "{sample}_R1_trimmed.fastq"), sample=SAMPLES),
        expand(pj(trim_trunc_path, "{sample}_R2_trimmed.fastq"), sample=SAMPLES),

        # Host-filtered
        expand(pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"), sample=SAMPLES),
        expand(pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz"), sample=SAMPLES),

        # Post-trim FastQC
        expand(pj(fastqc_path, "posttrim/{sample}", "{sample}_R1_trimmed.clean_1_fastqc.html"), sample=SAMPLES),
        expand(pj(fastqc_path, "posttrim/{sample}", "{sample}_R2_trimmed.clean_2_fastqc.html"), sample=SAMPLES),

        # MultiQC post-trim
        pj(fastqc_path, "posttrim_multiqc_report.html"),

        # Kraken2 + Bracken
        expand(pj(kraken_path, "{sample}.kraken.txt"), sample=SAMPLES),
        expand(pj(kraken_path, "{sample}.kreport2"), sample=SAMPLES),
        expand(pj(kraken_path, "{sample}.bracken.tsv"), sample=SAMPLES),

        # Aggregated taxonomy output
        pj(kraken_path, "Combined-taxonomy.tsv"),
        
        #Combined reads for HUMAnN
        expand(pj(nonhost_path, "{sample}_combined.fastq.gz"), sample=SAMPLES),
        
        #Genefamilies output from HUMAnN
        expand(pj(humann_path,"Individual","{sample}", "{sample}_genefamilies.tsv"), sample=SAMPLES),
        
        #Pathogen abundance estimates from HuMAnN
        expand(pj(humann_path,"Individual","{sample}", "{sample}_pathabundance.tsv"), sample=SAMPLES),

        #Pathogen coverage values from HuMAnN
        expand(pj(humann_path,"Individual","{sample}", "{sample}_pathcoverage.tsv"), sample=SAMPLES)
    
        
# ---------------------------------------------------------------------------------------------------------------------
#
# Processing and QC 
#
# ---------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------
# Unzip
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
# Pre-trim FASTQC
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
# MULTIQC Pre-trim
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
# Adapter Trimming
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
# Host Filtering
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
# Post-Trim FASTQC
# ------------------------------------------------------
rule fastqc_posttrim:
    input:
        r1=pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"),
        r2=pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz")
    output:
        dir=directory(pj(fastqc_path, "posttrim/{sample}")),
        r1_html=pj(fastqc_path, "posttrim/{sample}", "{sample}_R1_trimmed.clean_1_fastqc.html"),
        r2_html=pj(fastqc_path, "posttrim/{sample}", "{sample}_R2_trimmed.clean_2_fastqc.html")
    conda: "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {output.dir}
        fastqc {input.r1} {input.r2} -o {output.dir}
        """

# ------------------------------------------------------
# Post-trim MULTIQC
# ------------------------------------------------------
rule multiqc_posttrim:
    input:
        expand(pj(fastqc_path, "posttrim/{sample}"), sample=SAMPLES)
    output:
        html=pj(fastqc_path,"posttrim_multiqc_report.html")
    conda: "envs/fastqc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.html})
        multiqc {fastqc_path} -o $(dirname {output.html})
        mv $(dirname {output.html})/multiqc_report.html {output.html}
        """

# ---------------------------------------------------------------------------------------------------------------------
#
# Taxonomic Classification 
#
# ---------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------
# Kraken2 for read classification
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
    resources:
        mem_mb=300000,
        himem=True
    params:
        out_dir=kraken_path,
        database=kraken_db_loc
    threads: 16
    shell:
        """
        mkdir -p {params.out_dir}
        kraken2 \
            --no-mmap \
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
        reads=config["bracken_read_len"]
    threads: 4
    conda:
        "envs/kraken.yaml"
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
        combined=pj(kraken_path, "Combined-taxonomy.tsv")
    conda:
        "envs/kraken.yaml"
    shell:
        """
        echo -e "sample\tname\ttaxid\tabs\tfrac" > {output.combined}
        for f in {input}; do
            s=$(basename $f .bracken.tsv)
            awk -v s=$s 'NR>1{{print s "\\t" $0}}' $f >> {output.combined}
        done
        """

# ---------------------------------------------------------------------------------------------------------------------
#
# Functional Annotation 
#
# ---------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------
# Read Concatenation for HUMAnN
# ------------------------------------------------------
rule concat_nonhost_reads:
  input:
    r1_nonhost=pj(nonhost_path, "{sample}_R1_trimmed.clean_1.fastq.gz"),
    r2_nonhost=pj(nonhost_path, "{sample}_R2_trimmed.clean_2.fastq.gz")
  output:
    comb=pj(nonhost_path, "{sample}_combined.fastq.gz")
  resources:
    mem_mb=8000, # MB
    runtime=240 # 4 min
  shell:
    """
    cat {input.r1_nonhost} {input.r2_nonhost} > {output.comb}
    """
    
# ------------------------------------------------------
# HUMAnN Functional Annotation
# ------------------------------------------------------
rule run_humann_nonhost:
  input:
    # HUMAnN nucleotide database (ChocoPhlAn):
    # Used for species-resolved pangenome mapping of reads
    choc_db=config['chocoplhan_db_loc'],

    # HUMAnN protein database (UniRef):
    # Used for translated search of reads not mapped at the nucleotide level
    uniref_db=config['uniref_db_loc'],

    # MetaPhlAn Bowtie2 database:
    # Required for taxonomic profiling (MetaPhlAn 4, vJun23)
    metaphlan_db=config['metaphlan_db_loc'],
   
    # Concatenated, host-filtered paired-end reads (R1 + R2)
    # Input to HUMAnN for functional profiling
    comb_nonhost=pj(nonhost_path, "{sample}_combined.fastq.gz")
    
  output:
    # Gene family abundances (UniRef90 IDs)
    # Rows include both:
    #   - Stratified entries (gene|taxon)
    #   - Unstratified totals (gene only)
    genefams=pj(
        humann_path, "Individual", "{sample}",
        "{sample}_genefamilies.tsv"
    ),
    
    # Pathway abundance estimates
    # Quantifies the relative abundance of MetaCyc pathways
    # Includes stratified (taxon-specific) and unstratified values
    pathabund=pj(
        humann_path, "Individual", "{sample}",
        "{sample}_pathabundance.tsv"
    ),
    
    # Pathway coverage values
    # Represents pathway completeness (0–1), not abundance
    # Less sensitive to sequencing depth than pathabundance
    pathcov=pj(
        humann_path, "Individual", "{sample}",
        "{sample}_pathcoverage.tsv"
    ),

    # MetaPhlAn taxonomic profile used internally by HUMAnN
    # Lists detected microbial taxa and their relative abundances
    # Generated during the HUMAnN run and stored in the temp directory
    buglist=pj(
        humann_path, "Individual", "{sample}",
        "{sample}_humann_temp",
        "{sample}_metaphlan_bugs_list.tsv"
    )

  resources:
    mem_mb=64000,   # Memory allocation (64 GB)
    runtime=1500    # Runtime limit in minutes (~25 hours)

  threads: 32

  conda:
    "envs/humann.yaml"

  params:
    # Output directory for HUMAnN results
    out_dir=pj(humann_path, "Individual", "{sample}")

  shell:
    """
    mkdir -p {params.out_dir}

    humann -i {input.comb_nonhost} -o {params.out_dir} \
        --threads {threads} \
        --search-mode uniref90 \
        --nucleotide-database {input.choc_db} \
        --protein-database {input.uniref_db} \
        --metaphlan-options "--bowtie2db {input.metaphlan_db} \
                             --index mpa_vJan25_CHOCOPhlAnSGB_202503"
    """


# ------------------------------------------------------
# Agreggating Functional Data
# ------------------------------------------------------
rule aggregate_humann_outs:
    """
    Aggregate per-sample HUMAnN outputs into cohort-level tables and
    regroup UniRef90 gene families into biologically interpretable
    functional categories (MetaCyc, KEGG, GO, PFAM, EggNOG).
    """

    input:
        # Utility mapping database (used implicitly by HUMAnN regroup/rename)
        mapping_db=config['utility_mapping_db_loc'],

        # Force completion of all per-sample HUMAnN outputs before aggregation
        genefams=expand(
            pj(humann_path, "Individual", "{sample}", "{sample}_genefamilies.tsv"),
            sample=SAMPLES
        ),
        pathabund=expand(
            pj(humann_path, "Individual", "{sample}", "{sample}_pathabundance.tsv"),
            sample=SAMPLES
        ),
        pathcov=expand(
            pj(humann_path, "Individual", "{sample}", "{sample}_pathcoverage.tsv"),
            sample=SAMPLES
        ),

    output:
        # Joined taxonomic abundance tables
        pathabund_cons=pj(humann_path, "Consolidated", "all_pathabundance.tsv"),
        pathcov_cons=pj(humann_path, "Consolidated", "all_pathcoverage.tsv"),

        # Joined UniRef90 gene family table
        genefams_cons=pj(humann_path, "Consolidated", "all_genefamilies.tsv"),

        # MetaCyc reactions
        genefams_cons_rxn=pj(humann_path, "Consolidated", "all_genefamilies_rxn.tsv"),
        genefams_cons_rxn_named=pj(humann_path, "Consolidated", "all_genefamilies_rxn_named.tsv"),

        # KEGG orthologs, pathways, and modules
        genefams_cons_ko=pj(humann_path, "Consolidated", "all_genefamilies_ko.tsv"),
        genefams_cons_ko_named=pj(humann_path, "Consolidated", "all_genefamilies_ko_named.tsv"),
        genefams_cons_kp_named=pj(humann_path, "Consolidated", "all_genefamilies_kp_named.tsv"),
        genefams_cons_km_named=pj(humann_path, "Consolidated", "all_genefamilies_km_named.tsv"),

        # EggNOG
        genefams_cons_en=pj(humann_path, "Consolidated", "all_genefamilies_eggnog.tsv"),

        # Gene Ontology
        genefams_cons_go=pj(humann_path, "Consolidated", "all_genefamilies_go.tsv"),
        genefams_cons_go_named=pj(humann_path, "Consolidated", "all_genefamilies_go_named.tsv"),

        # PFAM
        genefams_cons_pfam=pj(humann_path, "Consolidated", "all_genefamilies_pfam.tsv"),
        genefams_cons_pfam_named=pj(humann_path, "Consolidated", "all_genefamilies_pfam_named.tsv"),

    resources:
        mem_mb=10000,   # ~10 GB RAM
        runtime=60,     # ~1 hour

    conda:
        "envs/humann.yaml"

    params:
        # Output directory for consolidated tables
        out_dir=pj(humann_path, "Consolidated"),

        # Directory containing per-sample HUMAnN results
        individual_dir=pj(humann_path, "Individual")

    shell:
        """
        # Create consolidated output directory
        mkdir -p {params.out_dir}

        ####################################################################
        # 1. Join per-sample HUMAnN outputs across all samples
        ####################################################################

        humann_join_tables \
            -i {params.individual_dir} \
            -o {output.pathabund_cons} \
            --file_name pathabundance.tsv \
            --search-subdirectories

        humann_join_tables \
            -i {params.individual_dir} \
            -o {output.pathcov_cons} \
            --file_name pathcoverage.tsv \
            --search-subdirectories

        humann_join_tables \
            -i {params.individual_dir} \
            -o {output.genefams_cons} \
            --file_name genefamilies.tsv \
            --search-subdirectories

        ####################################################################
        # 2. Regroup UniRef90 gene families into MetaCyc reactions
        ####################################################################

        humann_regroup_table \
            -i {output.genefams_cons} \
            -g uniref90_rxn \
            -o {output.genefams_cons_rxn}

        humann_rename_table \
            -i {output.genefams_cons_rxn} \
            -n metacyc-rxn \
            -o {output.genefams_cons_rxn_named}

        ####################################################################
        # 3. Regroup UniRef90 gene families into KEGG orthologs
        ####################################################################

        humann_regroup_table \
            -i {output.genefams_cons} \
            -g uniref90_ko \
            -o {output.genefams_cons_ko}

        humann_rename_table \
            -i {output.genefams_cons_ko} \
            -n kegg-orthology \
            -o {output.genefams_cons_ko_named}

        # KEGG pathways and modules are derived from KO-level tables
        humann_rename_table \
            -i {output.genefams_cons_ko} \
            -n kegg-pathway \
            -o {output.genefams_cons_kp_named}

        humann_rename_table \
            -i {output.genefams_cons_ko} \
            -n kegg-module \
            -o {output.genefams_cons_km_named}

        ####################################################################
        # 4. Regroup UniRef90 gene families into EggNOG
        ####################################################################

        humann_regroup_table \
            -i {output.genefams_cons} \
            -g uniref90_eggnog \
            -o {output.genefams_cons_en}

        ####################################################################
        # 5. Regroup UniRef90 gene families into Gene Ontology (GO)
        ####################################################################

        humann_regroup_table \
            -i {output.genefams_cons} \
            -g uniref90_go \
            -o {output.genefams_cons_go}

        humann_rename_table \
            -i {output.genefams_cons_go} \
            -n go \
            -o {output.genefams_cons_go_named}

        ####################################################################
        # 6. Regroup UniRef90 gene families into PFAM domains
        ####################################################################

        humann_regroup_table \
            -i {output.genefams_cons} \
            -g uniref90_pfam \
            -o {output.genefams_cons_pfam}

        humann_rename_table \
            -i {output.genefams_cons_pfam} \
            -n pfam \
            -o {output.genefams_cons_pfam_named}
        """
