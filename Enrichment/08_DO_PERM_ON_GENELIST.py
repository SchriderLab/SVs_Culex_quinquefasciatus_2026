#!/usr/bin/env python3
"""
Usage: python 08_DO_PERM_ON_GENELIST.py --real_sv real_seg_SVExtremeAF.tsv --permdir master_perms_AFExtreme  --ir_genes IR_genes.txt --n_perms 10000 --gene_buffer 1000
"""

import os
import sys
import argparse
import numpy as np


# Column positions in IR_genes.txt (0-indexed), after the header row
COL_CHROM = 0
COL_START = 1
COL_END = 2
COL_LOC_ID = 4   # "gene ID (loc)" column, e.g. LOC6032041
COL_DESC = 5     # "gene description" column


# -----------------------------
# Read a 4-column SV/segment file: chrom start end val
# -----------------------------
def read_segments(fname):
    rows = []
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, val = line.split("\t")[:4]
            rows.append((chrom, int(start), int(end), int(val)))
    return rows


# -----------------------------
# Load IR gene coordinates directly from IR_genes.txt
# Applies optional buffer to gene start/end coordinates.
# -----------------------------
def load_ir_genes(fname, gene_buffer=0):

    chrom_rows = {}
    descriptions = {}

    with open(fname) as f:
        first = True

        for line in f:
            line = line.rstrip("\n")

            if not line.strip():
                continue

            if first:
                first = False
                continue

            parts = line.split("\t")

            chrom = parts[COL_CHROM]

            # Apply buffer here
            start = int(parts[COL_START]) - gene_buffer
            end = int(parts[COL_END]) + gene_buffer

            gene_id = parts[COL_LOC_ID]
            desc = parts[COL_DESC] if len(parts) > COL_DESC else ""

            chrom_rows.setdefault(chrom, []).append(
                (start, end, gene_id)
            )

            descriptions[gene_id] = desc

    return chrom_rows, descriptions


# -----------------------------
# Build numpy index
# -----------------------------
def build_gene_index(chrom_rows):

    index = {}

    for chrom, rows in chrom_rows.items():

        starts = np.array([r[0] for r in rows])
        ends = np.array([r[1] for r in rows])
        names = np.array([r[2] for r in rows], dtype=object)

        index[chrom] = (starts, ends, names)

    return index


# -----------------------------
# Find genes overlapping SV segments
# -----------------------------
def genes_overlapping_sv(segments, gene_index):

    hits = set()

    for chrom, start, end, val in segments:

        if val != 1:
            continue

        gi = gene_index.get(chrom)

        if gi is None:
            continue

        g_start, g_end, g_name = gi

        mask = (start <= g_end) & (end >= g_start)

        if mask.any():
            hits.update(g_name[mask].tolist())

    return hits


# -----------------------------
# Read permutation SV files
# -----------------------------
def read_perm_segments(permdir, chroms, perm_id):

    rows = []

    for chrom in chroms:

        fpath = os.path.join(
            permdir,
            f"{chrom}_perm_{perm_id}.tsv"
        )

        if not os.path.exists(fpath):

            sys.stderr.write(
                f"WARNING: missing {fpath}, skipping this chrom "
                f"for perm {perm_id}\n"
            )

            continue

        rows.extend(read_segments(fpath))

    return rows


# -----------------------------
# Main
# -----------------------------
def main():

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    ap.add_argument(
        "--real_sv",
        required=True,
        help="genome_sv_segments.tsv from SV segmentation"
    )

    ap.add_argument(
        "--permdir",
        required=True,
        help="directory containing chromosome permutation files"
    )

    ap.add_argument(
        "--ir_genes",
        required=True,
        help="IR gene coordinate file"
    )

    ap.add_argument(
        "--n_perms",
        type=int,
        required=True
    )

    ap.add_argument(
        "--gene_buffer",
        type=int,
        default=0,
        help="bp added upstream and downstream of each gene "
             "before overlap testing (default: 0)"
    )

    ap.add_argument(
        "--chroms",
        default="NC_051861.1,NC_051862.1,NC_051863.1",
        help="comma-separated chromosome names"
    )

    ap.add_argument(
        "--out",
        default="ir_sv_enrichment_results.tsv",
        help="output file"
    )

    args = ap.parse_args()


    chroms = args.chroms.split(",")


    sys.stderr.write(
        "Loading insecticide resistance gene coordinates...\n"
    )

    sys.stderr.write(
        f"  Gene buffer: +/- {args.gene_buffer} bp\n"
    )


    chrom_rows, descriptions = load_ir_genes(
        args.ir_genes,
        gene_buffer=args.gene_buffer
    )


    total_ir = sum(len(v) for v in chrom_rows.values())

    sys.stderr.write(
        f"  {total_ir} IR genes loaded across "
        f"{len(chrom_rows)} chromosome(s)\n"
    )


    ir_index = build_gene_index(chrom_rows)


    # ---------------- REAL ----------------

    sys.stderr.write(
        "Reading real SV segmentation...\n"
    )

    real_segments = read_segments(args.real_sv)

    real_hits = genes_overlapping_sv(
        real_segments,
        ir_index
    )

    real_count = len(real_hits)


    print(f"Total IR genes tested: {total_ir}")

    print(
        f"Real IR genes overlapping an SV region: "
        f"{real_count} "
        f"({100 * real_count / total_ir:.1f}% of tested IR genes)"
    )

    print(
        "IR genes overlapping SV regions in real data:"
    )

    for g in sorted(real_hits):

        print(
            f"  {g}\t{descriptions.get(g, '')}"
        )


    # ---------------- PERMUTATIONS ----------------

    sys.stderr.write(
        "Running permutations...\n"
    )

    perm_counts = []


    for i in range(args.n_perms):

        perm_segments = read_perm_segments(
            args.permdir,
            chroms,
            i
        )

        perm_hits = genes_overlapping_sv(
            perm_segments,
            ir_index
        )

        perm_counts.append(
            len(perm_hits)
        )


        if (i + 1) % 50 == 0 or i == args.n_perms - 1:

            sys.stderr.write(
                f"  finished {i + 1}/{args.n_perms}\n"
            )


    perm_counts = np.array(perm_counts)

    mean_perm = perm_counts.mean()

    std_perm = perm_counts.std()


    n_ge = int(
        (perm_counts >= real_count).sum()
    )

    pval = n_ge / len(perm_counts)


    enrichment = (
        real_count / mean_perm
        if mean_perm != 0
        else 0
    )


    print("\n--- Permutation test summary ---")

    print(
        f"Permutations run: {len(perm_counts)}"
    )

    print(
        f"Mean IR genes overlapping SV in permuted data: "
        f"{mean_perm:.3f} (sd {std_perm:.3f})"
    )

    print(
        f"Real IR genes overlapping SV: {real_count}"
    )

    print(
        f"Fold enrichment (real / mean permuted): "
        f"{enrichment}"
    )

    print(
        f"Permutations with count >= real: "
        f"{n_ge}/{len(perm_counts)}"
    )

    print(
        f"P-value: {pval}"
    )


    with open(args.out, "w") as out:

        out.write(
            "perm_id\tir_genes_overlapping_sv\n"
        )

        for i, c in enumerate(perm_counts):

            out.write(
                f"{i}\t{c}\n"
            )

        out.write(
            f"REAL\t{real_count}\n"
        )


    sys.stderr.write(
        f"Per-permutation counts written to {args.out}\n"
    )


if __name__ == "__main__":
    main()
