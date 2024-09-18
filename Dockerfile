FROM snakemake/snakemake:stable
ADD environment.yaml .
RUN mamba create -n snakemake-tutorial --clone snakemake; \
    mamba env update -n snakemake-tutorial -f environment.yaml;
RUN mkdir -p /tmp/conda
ENV CONDA_PKGS_DIRS /tmp/conda

# FROM condaforge/mambaforge:latest
# LABEL io.github.snakemake.containerized="true"
# LABEL io.github.snakemake.conda_env_hash="b1d825247ab634fae2e4333c695a5d9736ef1783d7d5cb4f53de9346ca6945ec"

# # Step 1: Retrieve conda environments

# # Conda environment:
# #   source: envs/samtools.yaml
# #   prefix: /conda-envs/ce8b782c251483667307ba4ed3f68b5a
# #   channels:
# #     - bioconda
# #     - conda-forge
# #   dependencies:
# #     - samtools =1.9
# RUN mkdir -p /conda-envs/ce8b782c251483667307ba4ed3f68b5a
# COPY envs/samtools.yaml /conda-envs/ce8b782c251483667307ba4ed3f68b5a/environment.yaml

# # Step 2: Generate conda environments

# RUN mamba env create --prefix /conda-envs/ce8b782c251483667307ba4ed3f68b5a --file /conda-envs/ce8b782c251483667307ba4ed3f68b5a/environment.yaml && \
#     mamba clean --all -y
