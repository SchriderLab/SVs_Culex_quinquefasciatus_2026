#!/usr/bin/env python3
"""
Filter the master permuted SV files - do this for each perm so that only the SVs present in VCF of interest
are kept.
USAGE: ./03_filter_perms_by_dataset.py master_perms/ jasmine_overlapping_genes_high_snpeff_CLR_under100kb.vcf master_perms_filtbysnpclr
"""

import sys
import os
import glob
import re

CHROMS = {"NC_051861.1", "NC_051862.1", "NC_051863.1"}


def sv_key(chrom, pos, vcf_id, svtype, end):
    if vcf_id and vcf_id != ".":
        return vcf_id
    return f"{chrom}_{pos}_{end}_{svtype}"


def read_vcf_ids(fname):
    """read in and return every ID from the filtered real VCF"""
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
        print("Usage: filter_perms_by_dataset.py master_perm_dir dataset.vcf outdir")
        sys.exit(1)

    master_perm_dir, dataset_vcf, outdir = sys.argv[1:4]

    keep_ids = read_vcf_ids(dataset_vcf)
    print(f"{dataset_vcf}: {len(keep_ids)} real SVs to match against the master perms")

    os.makedirs(outdir, exist_ok=True)

    perm_files = glob.glob(os.path.join(master_perm_dir, "master_perm_*.tsv"))
    print(f"found {len(perm_files)} master permutation files")

    perm_re = re.compile(r"master_perm_(\d+)\.tsv$")

    matched_ever = set()

    for fpath in perm_files:
        m = perm_re.search(fpath)
        if not m:
            continue
        perm_idx = m.group(1)

        by_chrom = {}
        with open(fpath) as f:
            for line in f:
                sv_id, chrom, start, end = line.rstrip("\n").split("\t")
                if sv_id not in keep_ids:
                    continue
                matched_ever.add(sv_id)
                by_chrom.setdefault(chrom, []).append((start, end))

        for chrom, rows in by_chrom.items():
            outpath = os.path.join(outdir, f"{chrom}_perm_{perm_idx}.tsv")
            with open(outpath, "w") as out:
                for start, end in rows:
                    out.write(f"{chrom}\t{start}\t{end}\t1\n")

    unmatched = keep_ids - matched_ever
    print(f"done -- filtered permutations written to {outdir}")
    print(f"{len(matched_ever)}/{len(keep_ids)} real SV IDs found in the master set")
    if unmatched:
        print(f"WARNING: {len(unmatched)} SV(s) in {dataset_vcf} were not found "
              f"in the master permutations -- check that this dataset's VCF is "
              f"really a filtered subset of the master VCF (same IDs).")


if __name__ == "__main__":
    main()
