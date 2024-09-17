[![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/snakemake/snakemake-tutorial-data)

# Example data for the Snakemake tutorial

This repository hosts the data needed for the [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

## Notes
* Rules are tasks to be carried out and must have an input, output, and script

Snakemake automatically runs a file named Snakefile but you can run other files using the -s option:
```
snakemake -s <file> --cores 1
```
