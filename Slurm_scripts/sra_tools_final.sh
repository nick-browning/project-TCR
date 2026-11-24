#!/bin/bash
#SBATCH --job-name=sra_tools_final.sh
#SBATCH --partition=bigmem
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

OUTDIR=/data/scratch/projects/punim1612/nick/sra_tools_sample_run_2_BCR
ACC_LIST=/data/projects/punim1612/CSI/Nick-project/Slurm_scripts/sample_accessions_BCR.txt

#pysradb 2.5.1 to grab metadata for SRA runs, then can rename R1/R2 based on data type for correct input into cellrangermulti_config.csv

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

    # First pass: measure the modal read length of every split file. Index
    # reads (<20 bp) are removed; the rest are kept as candidates in file order.
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

    # Second pass: assign R1/R2 relative to each other rather than by fixed
    # cutoffs. Some runs are sequenced 150x150, so R1 (barcode+UMI+read-through)
    # is the same length as R2. The shorter read is R1; on a tie, fasterq-dump
    # file order (_1 before _2) decides.
    if (( ${#reads[@]} != 2 )); then
        echo "ERROR: expected 2 non-index reads for $SRR, found ${#reads[@]} - leaving files unrenamed" >&2
        continue
    fi
    if (( lens[1] < lens[0] )); then
        r1=${reads[1]}; r2=${reads[0]}
    else
        r1=${reads[0]}; r2=${reads[1]}
    fi
    echo "  - calling R1: $(basename "$r1"), R2: $(basename "$r2")"
    # -n: never overwrite an existing file, had an error in a previous run where R1 was overridden by R2
    mv -n "$r1" ${OUTDIR}/${SRR}_S1_L001_R1_001.fastq.gz
    mv -n "$r2" ${OUTDIR}/${SRR}_S1_L001_R2_001.fastq.gz

done < "$ACC_LIST"
