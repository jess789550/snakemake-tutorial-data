configfile: "config.yaml"

# SAMPLES = ["A", "B"]

# target rule with default output target for command line
rule all:
    input:
        "plots/quals.svg"

# snakemake -n

# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/A.fastq"
#     output:
#         "mapped_reads/A.bam"
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"

# snakemake -np mapped_reads/A.bam  # -n is dry run and -p is print command
# snakemake --cores 1 mapped_reads/A.bam  # need to give target output on command line and number of cores

# rule bwa_map:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq" # generalise the FASTQ file
#     output:
#         "mapped_reads/{sample}.bam" # generalise the BAM file
#     shell:
#         "bwa mem {input} | samtools view -Sb - > {output}"

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

# rule bwa_map:
#     input:
#         "data/genome.fa",
#         get_bwa_map_input_fastqs
#     output:
#         "mapped_reads/{sample}.bam"
#     threads: 8
#     shell:
#         "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

# rule bwa_map:
#     input:
#         "data/genome.fa",
#         get_bwa_map_input_fastqs
#     output:
#         "mapped_reads/{sample}.bam"
#     params:
#         rg=r"@RG\tID:{sample}\tSM:{sample}"
#     threads: 8
#     shell:
#         "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"

rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

# snakemake -np mapped_reads/B.bam
# snakemake -np mapped_reads/A.bam mapped_reads/B.bam
# snakemake -np mapped_reads/{A,B}.bam

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

# snakemake -np sorted_reads/B.bam

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

# snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg

# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#     output:
#         "calls/all.vcf"
#     shell:
#         "bcftools mpileup -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        mutation_rate=config["prior_mutation_rate"] 
    log:
        "logs/bcftools_call.log" # Log files must contain exactly the same wildcards as the output files to avoid file name clashes between different jobs of the same rule.
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -P {params.mutation_rate} -mv - > {output} 2> {log}"

# snakemake --dag calls/all.vcf | dot -Tsvg > dag2.svg
# snakemake -n --forcerun bcftools_call

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"

# snakemake --dag calls/all.vcf | dot -Tsvg > dag3.svg
# snakemake --cores 1 plots/quals.svg
# snakemake --help
# snakemake --forcerun samtools_sort --cores 1
# snakemake --forcerun samtools_index --cores 1 -np
# snakemake --cores 1 --summary # don't use --summary because there's a bug
