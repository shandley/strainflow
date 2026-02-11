#!/usr/bin/env bash
set -euo pipefail

REF_DIR="${1:-data/reference}"

echo "Creating reference directories..."
mkdir -p "${REF_DIR}/human_genome"
mkdir -p "${REF_DIR}/species_genomes"

# Download human genome for host removal
echo "Downloading GRCh38 human genome..."
if [ ! -f "${REF_DIR}/human_genome/GRCh38.fa.gz" ]; then
  wget -O "${REF_DIR}/human_genome/GRCh38.fa.gz" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
fi

echo "Building bowtie2 index for human genome..."
if [ ! -f "${REF_DIR}/human_genome/GRCh38.1.bt2" ]; then
  gunzip -k "${REF_DIR}/human_genome/GRCh38.fa.gz"
  bowtie2-build --threads 8 "${REF_DIR}/human_genome/GRCh38.fa" "${REF_DIR}/human_genome/GRCh38"
fi

echo "Reference setup complete."
echo ""
echo "NOTE: You still need to provide species reference genomes."
echo "Place concatenated reference FASTA at: ${REF_DIR}/species_genomes/all_species.fasta"
echo "Build bowtie2 index: bowtie2-build all_species.fasta all_species_bt2"
echo "Create scaffold-to-bin file: ${REF_DIR}/species_genomes/scaffold_to_bin.stb"
