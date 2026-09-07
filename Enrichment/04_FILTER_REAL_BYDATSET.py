#!/usr/bin/env python3

"""
Get the real SV interval file (not permuted) from the master SV list.
USAGE:python 04_FILTER_REAL_BYDATSET.py master_svs.tsv jasmine_overlapping_genes_high_snpeff_CLR_under100kb.vcf out_real_segments.tsv
"""

import sys

CHROMS = {"NC_051861.1", "NC_051862.1", "NC_051863.1"}


def sv_key(chrom, pos, vcf_id, svtype, end):
    if vcf_id and vcf_id != ".":
        return vcf_id
    return f"{chrom}_{pos}_{end}_{svtype}"


def read_vcf_ids(fname):
    """SV IDs present in this dataset's own real (filtered) VCF."""
    ids = set()
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, vcf_id, info = fields[0], fields[1], fields[2], fields[7]

            if chrom not in CHROMS:
                continue

            svtype = ""
            end = ""
            for entry in info.split(";"):
                if entry.startswith("SVTYPE="):
                    svtype = entry.split("=", 1)[1]
                elif entry.startswith("END="):
                    end = entry.split("=", 1)[1]

            # insertions: no reference span, so treat as a single-bp point at
            # POS (same convention as 01_extract_master_svs.sh)
            if svtype == "INS":
                end = pos
            else:
                if end == "" or end == pos:
                    continue

            ids.add(sv_key(chrom, pos, vcf_id, svtype, end))

    return ids


def main():
    if len(sys.argv) != 4:
        print("Usage: filter_real_by_dataset.py master_svs.tsv dataset.vcf "
              "out_real_segments.tsv")
        sys.exit(1)

    master_file, dataset_vcf, out_file = sys.argv[1:4]

    keep_ids = read_vcf_ids(dataset_vcf)
    print(f"{dataset_vcf}: {len(keep_ids)} real SVs to match against the master list")

    n_written = 0
    matched_ids = set()

    with open(master_file) as f, open(out_file, "w") as out:
        for line in f:
            sv_id, chrom, start, end = line.rstrip("\n").split("\t")
            if sv_id not in keep_ids:
                continue
            matched_ids.add(sv_id)
            out.write(f"{chrom}\t{start}\t{end}\t1\n")
            n_written += 1

    unmatched = keep_ids - matched_ids
    print(f"wrote {n_written} real SV intervals to {out_file}")
    if unmatched:
        print(f"WARNING: {len(unmatched)} SV(s) in {dataset_vcf} were not found "
              f"in the master list -- check that this dataset's VCF is really "
              f"a filtered subset of the master VCF (same IDs).")


if __name__ == "__main__":
    main()
