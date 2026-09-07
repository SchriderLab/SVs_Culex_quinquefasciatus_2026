#!/usr/bin/env python3

"""
Shift the master sv file once per permutation and make sure interval stays in tact

A circular shift that uses a single random offset per chrom pre permutations gets applied to every interval on the chrom
This preserves the relative spacing between the SVs

USAGE: ./02_PERMUTE_MASTER_SV.py master_svs.tsv GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna.fai master_perms 10000
"""

import sys
import os
import random

"""read in that genome.fai file (extract the size of chrom because you will rotate the chrom based on a random number from 0-size ofchrom)"""
def read_fai(fname):
    lengths = {}
    with open(fname) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom, length = parts[0], int(parts[1])
            lengths[chrom] = length
    return lengths

"""read in the master SV file"""
def read_master(fname):
    by_chrom = {}
    with open(fname) as f:
        for line in f:
            if not line.strip():
                continue
            sv_id, chrom, start, end = line.rstrip("\n").split("\t")
            by_chrom.setdefault(chrom, []).append((sv_id, int(start), int(end)))
    return by_chrom

"""the shifting function, if sv end is less than the end of chrom, just print, otherwise wrap and split it and prerserve the OG ID for both intervals"""
def shift_svs(svs, chrom_len):

    shift = random.randint(0, chrom_len - 1)
    shifted = []

    for sv_id, start, end in svs:
        length = end - start
        new_start = ((start - 1 + shift) % chrom_len) + 1
        new_end = new_start + length

        if new_end <= chrom_len:
            shifted.append((sv_id, new_start, new_end))
        else:
            shifted.append((sv_id, new_start, chrom_len))
            shifted.append((sv_id, 1, new_end - chrom_len))

    return shifted


def main():
    if len(sys.argv) != 5:
        print("Usage: permute_master_svs.py master_svs.tsv genome.fai outdir n_perms")
        sys.exit(1)

    master_file, fai_file, outdir, n_perms = sys.argv[1:5]
    n_perms = int(n_perms)

    chrom_len = read_fai(fai_file)
    by_chrom = read_master(master_file)

    os.makedirs(outdir, exist_ok=True)

    print("chromosomes:", list(by_chrom.keys()))
    print("permutations:", n_perms)

    for p in range(n_perms):
        outpath = os.path.join(outdir, f"master_perm_{p}.tsv")
        with open(outpath, "w") as out:
            for chrom, svs in by_chrom.items():
                if chrom not in chrom_len:
                    raise ValueError(f"{chrom} missing from .fai file")
                shifted = shift_svs(svs, chrom_len[chrom])
                for sv_id, s, e in shifted:
                    out.write(f"{sv_id}\t{chrom}\t{s}\t{e}\n")

        if (p + 1) % 10 == 0 or p == n_perms - 1:
            print(f"finished {p + 1}/{n_perms}")

    print("done")


if __name__ == "__main__":
    main()
