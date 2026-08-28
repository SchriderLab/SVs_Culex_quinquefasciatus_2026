#1. we followed the 0. Quality Control and Read Mapping and 1. Variant Calling and Filtering workflows from Haba et al. (2025), found in https://github.com/YukiHaba/PipPop_molestus_origin/tree/main/scripts
# with the following modifications:
  # --include 'F_MISSING < 0.25 & MQ > 40 & QUAL > 30' \ in 1. variant_calling.sh was changed to --include 'F_MISSING < 0.5 & MQ > 40 & QUAL > 30' \
  # we ran the scripts only on the samples identified as C. quinquefasciatus in the original study

#2. subset VCFs by population (South Africa, Madagascar, Cameroon/Gabon)
bcftools view -S safrica.txt good_biallelic_snps_05.rm.combined.accessible.vcf.gz -Oz -o snps_final_southafrica.vcf.gz

#3. run vcf2sf.py (https://github.com/SchriderLab/timesweeper-experiments/blob/main/scripts/comp_methods/vcf2sf.py) - generates the "FreqFile" for each chromosome to use for SweepFinder

#4. SweepFinder commands:
# "CombinedFreqFile" = the 3 "FreqFiles" for each chromosome in a population concatenated together

# run the following for each population to generate SFS:
SweepFinder2 -f CombinedFreqFile SpectFile

# run the following for each population. "OutFile" = CLR results
SweepFinder2 -lg 1000 FreqFile SpectFile OutFile

# we parallelized by splitting the FreqFiles into overlapping chunks, e.g.
CHUNK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chunks.txt)
SweepFinder2 -lg 1000 southafrica_chunks/chr1/${CHUNK} SAfrica_SpectFile.txt southafrica_chunks/chr1/sweepfinder_southafrica_chr1_${CHUNK}

# to generate these chunks:
OUTPREFIX="chunk"
NCHUNKS=10
OVERLAP_BP=50000   # overlap in basepairs

# Get min and max genomic positions
MIN_POS=$(awk 'NR==2 {print $1}' "$INPUT")
MAX_POS=$(awk 'END {print $1}' "$INPUT")

RANGE=$((MAX_POS - MIN_POS))
WINDOW_SIZE=$((RANGE / NCHUNKS))

echo "Min position: $MIN_POS"
echo "Max position: $MAX_POS"
echo "Window size: $WINDOW_SIZE"
echo "Overlap (bp): $OVERLAP_BP"

HEADER=$(head -n 1 "$INPUT")

for ((i=0; i<NCHUNKS; i++)); do

    START=$((MIN_POS + i * WINDOW_SIZE))
    END=$((MIN_POS + (i + 1) * WINDOW_SIZE))

    # Add overlap (except first and last)
    if [ $i -ne 0 ]; then
        START=$((START - OVERLAP_BP))
    fi

    if [ $i -ne $((NCHUNKS - 1)) ]; then
        END=$((END + OVERLAP_BP))
    else
        END=$MAX_POS
    fi

    OUTFILE="chr3/${OUTPREFIX}_$((i+1)).txt"

    echo "Writing $OUTFILE (positions $START to $END)"

    {
        echo "$HEADER"
        awk -v s="$START" -v e="$END" 'NR>1 && $1 >= s && $1 <= e' "$INPUT"
    } > "$OUTFILE"

done

# then concatenate results:
# Keep header from first chunk
head -n 1 chr3/sweepfinder_madagascar_chr3_chunk_1.txt > sweepfinder_madagascar_chr3_merged.txt

# Concatenate all data rows, sort numerically by location, remove duplicates
tail -n +2 chr3/sweepfinder_madagascar_chr3_chunk_*.txt \
    | sort -n -k1,1 \
    | awk '!seen[$1]++' \
    >> sweepfinder_madagascar_chr3_merged.txt

