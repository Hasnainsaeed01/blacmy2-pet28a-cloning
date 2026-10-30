#!/usr/bin/env bash
# Re-fetch the input sequences for the blaCMY-2 / pET-28a(+) cloning project.
# The sequences are also bundled in data/ (they are small); this script shows
# exactly where each one came from so the project is fully reproducible.
set -euo pipefail
cd "$(dirname "$0")"

echo "Downloading pET-28a(+) vector (5369 bp, NCBI PU845465)..."
curl -s --max-time 60 \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=3176436104&rettype=fasta&retmode=text" \
  -o pET28a_vector_downloaded.fasta
echo "Vector saved to pET28a_vector_downloaded.fasta"

echo ""
echo "blaCMY-2 comes from the user's own Salmonella Newport SL254 assembly"
echo "(GCF_000016045.1, plasmid NC_009140.1, plus strand 72367-73512,"
echo "complete 1146 bp ORF). To re-derive it from scratch:"
echo "  1. datasets download genome accession GCF_000016045.1 --include genome"
echo "  2. Extract NC_009140.1:72367-73512 (+) and verify ATG start / TAA stop:"
echo "     python3 -c \"from Bio import SeqIO; ...\""
echo ""
echo "The bundled data/blaCMY2_insert.fasta is that exact 1146 bp sequence."
