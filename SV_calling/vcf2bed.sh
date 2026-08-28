awk '
BEGIN { OFS="\t" }
!/^#/ {
    chrom = $1
    start = $2 - 1
    id = $3

    end = start
    if (match($8, /END=([0-9]+)/, arr)) {
        end = arr[1]
    }

    print chrom, start, end, id
}
' in.vcf  > out.bed
