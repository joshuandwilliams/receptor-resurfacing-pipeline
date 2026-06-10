#!/bin/bash
#SBATCH --job-name="AF3_Pwt3"
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem=64G
#SBATCH -p jic-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=AF3_Pwt3.out
#SBATCH --error=AF3_Pwt3.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

# Input: PWT3 mature chain (residues 19-141 of UniProt A0A223ZP76, 123 aa).
# A0A223ZP76 is the B-type (virulent, Br48 isolate). The Ao-type avirulent allele
# (DDBJ LC202650 etc.) is the actual design target per the lit sweep; the two differ
# at 12 aa, which does not change the fold but matters for eventual binder design.
# Swap the JSON sequence for the Ao-type when re-running for design.

export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_FORCE_UNIFIED_MEMORY=true
export XLA_CLIENT_MEM_FRACTION=3.2

AF3_MODEL_DIR="/hpc-home/jowillia/singularity/AlphaFold3"
AF3_DATA_DIR="${HOME}/singularity/AlphaFold3/af3_db"
JSON_INPUT="data/pwt3_monomer_af3.json"
OUTPUT_DIR="./alphafold3_pwt3_output"

mkdir -p "${OUTPUT_DIR}"

if [ ! -f "${JSON_INPUT}" ]; then
    exit 1
fi

SEQ_LEN=$(python3 -c "import json; d=json.load(open('${JSON_INPUT}')); print(len(d['sequences'][0]['protein']['sequence']))")

mkdir -p "${AF3_DATA_DIR}"
for f in /nbi/Reference-Data/AlphaFold/db-v3.0.0/*; do
    ln -sfn "$f" "${AF3_DATA_DIR}/$(basename "$f")"
done

ln -sfn /nbi/Reference-Data/AlphaFold/db-v2.3.2/small_bfd/bfd-first_non_consensus_sequences.fasta \
    "${AF3_DATA_DIR}/bfd-first_non_consensus_sequences.fasta"
ln -sfn /nbi/Reference-Data/AlphaFold/db-v2.3.2/mgnify/mgy_clusters_2022_05.fa \
    "${AF3_DATA_DIR}/mgy_clusters_2022_05.fa"

source package e8edb411-7374-4342-b9f1-408da41fc197

srun run_alphafold.py \
    --json_path="${JSON_INPUT}" \
    --model_dir="${AF3_MODEL_DIR}" \
    --db_dir="${AF3_DATA_DIR}" \
    --output_dir="${OUTPUT_DIR}"
