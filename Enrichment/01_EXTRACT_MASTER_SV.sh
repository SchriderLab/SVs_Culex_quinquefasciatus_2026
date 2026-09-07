#!/bin/bash
set -euo pipefail

#Makes a SV interval list that includes id,chrom,start,end
#Keeps each SV as its own row, keeps ID for downstream filtering
#USAGE: ./01_EXTRACT_MASTER_SV.sh jasmine_filt_no0_under100kb.vcf master_svs.tsv

VCF="$1"
OUT="$2"

grep -v "^#" "$VCF" | awk -F'\t' '
{
    chr = $1
    pos = $2
    id  = $3
    info = $8

    if (chr != "NC_051861.1" && chr != "NC_051862.1" && chr != "NC_051863.1")
        next

    svtype = ""
    end = ""

    n = split(info, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], b, "=")
        if (b[1] == "SVTYPE") svtype = b[2]
        if (b[1] == "END")    end = b[2]
    }

    # insertions: no reference span, so treat as a single-bp point at POS.
    # (END is usually already == POS per the VCF spec, but set it explicitly
    # in case END is missing/malformed for some INS records)
    if (svtype == "INS") {
        end = pos
    } else {
        # non-insertion SVs still need a real interval
        if (end == "" || end == pos) next
    }

    key = (id == "" || id == ".") ? (chr "_" pos "_" end "_" svtype) : id

    print key "\t" chr "\t" pos "\t" end
}
' | sort -k1,1 > "$OUT"

n=$(wc -l < "$OUT")
echo "Wrote $n master SV intervals to $OUT"
