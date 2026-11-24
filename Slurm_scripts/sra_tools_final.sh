git #!/bin/bash
#SBATCH --job-name=sra_tools_final.sh
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G                     
#SBATCH --time=24:00:00               
#SBATCH --output=sra_tools%j.out
#SBATCH --error=sra_tools%j.err

# Download SRA runs and convert to paired-end FASTQ
# accessions.txt: one SRA/ERR/DRR run accession per line (blank lines and
# lines starting with # are skipped).

# Load prerequisite modules, then load SRA-tools.
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load SRA-Toolkit/3.0.5

#Hardcoded directories and files for this script, will need to change to nextflow parameters in future
OUTDIR=/data/scratch/projects/punim1612/nick/sra_tools_sample_run_2_BCR
ACC_LIST=/data/projects/punim1612/CSI/Nick-project/Slurm_scripts/sample_accessions_BCR.txt
SCRIPTS=/data/projects/punim1612/CSI/Nick-project/Slurm_scripts
DETECT_R1=$SCRIPTS/detect_r1.py
WL_DIR=/apps/easybuild-2022/easybuild/software/Core/CellRanger/9.0.1/lib/python/cellranger/barcodes

#pysradb 2.5.1 to grab metadata for SRA runs, then can rename R1/R2 based whitelist match for correct input into cellrangermulti_config.csv

for SRR in $(cat "$ACC_LIST"); do
    echo "Processing $SRR..."

    TMPDIR="$OUTDIR/temp/${SRR}_temp"
    mkdir -p "$TMPDIR"

    prefetch "$SRR" -O "$TMPDIR" --max-size 1000g </dev/null

    # Point fasterq-dump at the prefetched run dir ($TMPDIR/$SRR/$SRR.sra),
    # otherwise it ignores the local copy and re-streams from NCBI.
    fasterq-dump "$TMPDIR/$SRR" --include-technical --split-files --outdir "$TMPDIR" --temp "$TMPDIR" --threads 8 </dev/null

    gzip "$TMPDIR"/*.fastq
    mv "$TMPDIR"/*.fastq.gz "$OUTDIR"/
    rm -rf "$TMPDIR"

    # First pass: measure the modal read length of every split file. Index / technical
    # reads (<20 bp) are removed; the rest are kept as candidates in file order, passed into detect_r1.py for 
    # assignemnt of read 1 vs read 2.
    reads=()
    lens=()
    for i in $(seq 1 4); do
        f=${OUTDIR}/${SRR}_${i}.fastq.gz
        [[ -f "$f" ]] || continue
        temp_len=$(zcat "$f" | head -n400 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | sort -rn | head -n 1 | awk '{print $2}')
        echo "${SRR}_${i}: read length is ${temp_len}"
        if (( temp_len >= 20 )); then
            reads+=("$f")
            lens+=("$temp_len")
        else
            echo "  - index read, removing"
            rm "$f"
        fi
    done

    if (( ${#reads[@]} != 2 )); then
        echo "ERROR: expected 2 non-index reads for $SRR, found ${#reads[@]} - leaving files unrenamed" >&2
        continue
    fi



done < "$ACC_LIST"
