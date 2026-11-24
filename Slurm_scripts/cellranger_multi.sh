#!/bin/bash
#SBATCH --job-name=cellranger_multi
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=cellranger_multi%j.out
#SBATCH --error=cellranger_multi%j.err

# Run cellranger multi on a config CSV written by ../cellranger_csv_generator.py.
# Any mix of GEX / TCR / CITE / FB / BEAM libraries works; the generator only
# writes the [gene-expression] / [vdj] / [feature] sections the libraries need.
#
# Usage: sbatch cellranger_multi.sh [config.csv] [run_id]

set -euo pipefail

CSV=${1:-/data/projects/punim1612/CSI/Nick-project/Slurm_scripts/cellrangermulti_config.csv}
RUN_ID=${2:-multi_$(date +%Y%m%d_%H%M)}
OUTROOT=/data/scratch/projects/punim1612/nick/Project/cellranger_multi_output_sample_run_2

# Refuse to run an unfilled template (e.g. <GEX_SAMPLE_NAME> placeholders)
if grep -q '<' "$CSV"; then
    echo "ERROR: $CSV still has placeholders - write it with cellranger_csv_generator.py" >&2
    exit 1
fi

# cellranger refuses to run into an existing output directory
if [[ -e "$OUTROOT/$RUN_ID" ]]; then
    echo "ERROR: $OUTROOT/$RUN_ID already exists - pick a new run_id" >&2
    exit 1
fi

mkdir -p "$OUTROOT"

module load CellRanger/9.0.1

cellranger multi --id="$RUN_ID" --csv="$CSV" \
    --output-dir="$OUTROOT/$RUN_ID" \
    --localcores="${SLURM_CPUS_PER_TASK:-16}" --localmem=120
