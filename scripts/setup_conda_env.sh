#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-strainflow_env}"

echo "Creating conda environment: ${ENV_NAME}"
conda create -n "${ENV_NAME}" python=3.10 -y

echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

echo "Installing bioinformatics tools..."
conda install -c bioconda -c conda-forge \
  bowtie2 \
  samtools \
  fastp \
  megahit \
  blast \
  -y

echo "Installing Python packages..."
pip install inStrain metaphlan

echo "Verifying installations..."
fastp --version
bowtie2 --version | head -1
samtools --version | head -1
megahit --version
blastn -version | head -1
inStrain --version
metaphlan --version

echo "Environment ${ENV_NAME} setup complete."
