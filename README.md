[![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/snakemake/snakemake-tutorial-data)

# Example data for the Snakemake tutorial

This repository hosts the data needed for the [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

## Notes
* Rules are tasks to be carried out and must have an input, output, and a shell script

Snakemake automatically runs a file named Snakefile but you can run other files using the -s option:
```
snakemake -s <file> --cores 1
```

Alternatively you can use modularisation:
```
include: "path/to/other.smk"
```

Snakemake wrappers repository: 
https://github.com/snakemake/snakemake-wrappers

Snakemake can be containerised in Docker:
```
snakemake --containerize > Dockerfile
```

Also jobs can be run in containers:
```
rule NAME:
    input:
        "table.txt"
    output:
        "plots/myplot.pdf"
    container:
        "docker://joseespinosa/docker-r-ggplot2:1.0"
    script:
        "scripts/plot-stuff.R"
```
