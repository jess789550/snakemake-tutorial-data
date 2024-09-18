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

# rule bwa_map:
#     input:
#         "data/genome.fa",
#         get_bwa_map_input_fastqs
#     output:
#         "mapped_reads/{sample}.bam"
#     params:
#         rg=r"@RG\tID:{sample}\tSM:{sample}"
#     log:
#         "logs/bwa_mem/{sample}.log"
#     threads: 8
#     shell:
#         "(bwa mem -R '{params.rg}' -t {threads} {input} | "
#         "samtools view -Sb - > {output}) 2> {log}"

rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        temp("mapped_reads/{sample}.bam") # deleted after pipeline is finished
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt" # duration of job
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"


# rule bwa_mem:
#   input:
#       ref="data/genome.fa",
#       sample=lambda wildcards: config["samples"][wildcards.sample]
#   output:
#       temp("mapped_reads/{sample}.bam")
#   log:
#       "logs/bwa_mem/{sample}.log"
#   params:
#       "-R '@RG\tID:{sample}\tSM:{sample}'"
#   threads: 8
#   wrapper:
#       "0.15.3/bio/bwa/mem" # A wrapper is a short script that wraps (typically) a command line application and makes it directly addressable from within Snakemake. For this, Snakemake provides the wrapper directive that can be used instead of shell, script, or run.


# snakemake -np mapped_reads/B.bam
# snakemake -np mapped_reads/A.bam mapped_reads/B.bam
# snakemake -np mapped_reads/{A,B}.bam
# snakemake --cores 1 mapped_reads/A.bam # target files are not deleted

# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam") # cannot be overwritten
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

# snakemake -np sorted_reads/B.bam

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
      "envs/samtools.yaml"
    shell:
        "samtools index {input}"

# snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg
# snakemake --software-deployment-method conda --cores 1
# snakemake --sdm conda -c 1

