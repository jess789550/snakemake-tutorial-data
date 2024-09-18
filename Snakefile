configfile: "config.yaml"

# SAMPLES = ["A", "B"]

# target rule with default output target for command line
rule all:
    input:
        "plots/quals.svg"

# snakemake -n

include: "Readmapping.smk"

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
        "logs/bcftools_call/all.log" # Log files must contain exactly the same wildcards as the output files to avoid file name clashes between different jobs of the same rule.
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
# snakemake --cores 1 --forceall
