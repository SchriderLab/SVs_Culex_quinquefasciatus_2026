#!/usr/bin/env python3
"""
GO enrichment via permutation testing on SV intervals.
Usage:python3 05_GET_TERM_COUNTS.py --real real_seg_outlier.tsv --permdir master_perms_outlier --gff GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.gff --go cleaned_go.tsv --n_perms 10000 --term_outdir permutations_out_outlier
"""

import os
import argparse
import numpy as np
import pandas as pd
from collections import defaultdict


# -----------------------------
# Load SV intervals
# -----------------------------
def read_intervals(fname):
    svs = []
    with open(fname) as f:
        for line in f:
            if not line.strip():
                continue
            chrom, start, end, val = line.strip().split()[:4]
            chrom = chrom.strip()
            svs.append((chrom, int(start), int(end), float(val)))
    return svs


# -----------------------------
# Load genes from GFF
# -----------------------------
def load_genes(gff_file):
    gff = pd.read_csv(gff_file, sep="\t", comment="#", header=None)

    gff.columns = [
        "seqid", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ]

    gff = gff[gff["type"] == "gene"].copy()

    def get_gene(attr):
        if "gene=" in attr:
            return attr.split("gene=")[1].split(";")[0]
        if "ID=gene-" in attr:
            return attr.split("ID=gene-")[1].split(";")[0]
        return None

    gff["gene"] = gff["attributes"].apply(get_gene)

    return gff[["seqid", "start", "end", "gene"]]


# -----------------------------
# Precomputed per-chromosome gene index (numpy arrays) for fast overlap
# -----------------------------
def build_gene_index(genes):
    index = {}
    for chrom, sub in genes.groupby("seqid"):
        index[chrom] = {
            "start": sub["start"].to_numpy(),
            "end": sub["end"].to_numpy(),
            "name": sub["gene"].to_numpy(dtype=object),
        }
    return index


# -----------------------------
# Interval overlap (vectorized)
# -----------------------------
def overlap(svs, gene_index):
    """
    svs: list of (chrom, start, end, val)
    gene_index: output of build_gene_index() -- built ONCE and reused across
                the real data and every permutation.
    Returns: list of (chrom, start, end, val, hit_gene_set), same order/filtering
             semantics as the original (SVs with val <= 0 are dropped).
    """
    by_chrom = defaultdict(list)
    for idx, (chrom, s, e, val) in enumerate(svs):
        if val <= 0:
            continue
        by_chrom[chrom].append((idx, s, e, val))

    results = [None] * len(svs)

    for chrom, rows in by_chrom.items():
        gi = gene_index.get(chrom)

        if gi is None or len(gi["start"]) == 0:
            for idx, s, e, val in rows:
                results[idx] = (chrom, s, e, val, set())
            continue

        g_start = gi["start"]
        g_end = gi["end"]
        g_name = gi["name"]

        sv_idx = np.array([r[0] for r in rows])
        sv_s = np.array([r[1] for r in rows])
        sv_e = np.array([r[2] for r in rows])
        sv_val = [r[3] for r in rows]

        # (n_sv, n_genes) boolean overlap matrix, computed in one vectorized shot
        mask = (sv_s[:, None] <= g_end[None, :]) & (sv_e[:, None] >= g_start[None, :])

        for row_i in range(len(rows)):
            hit_names = g_name[mask[row_i]]
            hit_set = {n for n in hit_names if n is not None}
            results[sv_idx[row_i]] = (chrom, sv_s[row_i], sv_e[row_i], sv_val[row_i], hit_set)

    return [r for r in results if r is not None]


# -----------------------------
# Read GO annotations
# -----------------------------
def read_go(fname):
    gene2go = {}
    with open(fname) as f:
        for line in f:
            if not line.strip():
                continue
            gene, go, _ = line.strip().split("\t")[:3]
            gene2go.setdefault(gene, set()).add(go)
    return gene2go


def build_tested_go_cache(genes, gene2go):
    """Compute the set of GO terms that are 'tested' (annotated on at least one
    gene in the GFF) ONCE, instead of recomputing it every permutation."""
    tested = set()
    for g in genes["gene"]:
        if g in gene2go:
            tested.update(gene2go[g])
    return tested


def fresh_tested_go_sets(tested_go_terms):
    """Cheap per-iteration reset: build a fresh {go: set()} dict from the cached
    term list instead of re-walking all genes/gene2go every time."""
    return {go: set() for go in tested_go_terms}


# -----------------------------
# Optional GO term name lookup
# -----------------------------
def read_go_names(fname):
    names = {}
    with open(fname) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                names[parts[0]] = parts[1]
    return names


# -----------------------------
# Convert SV hits -> GO gene sets
# -----------------------------
def go_gene_sets(sv_hits, gene2go):
    go2genes = {}
    for _, _, _, _, genes in sv_hits:
        for g in genes:
            if g not in gene2go:
                continue
            for go in gene2go[g]:
                go2genes.setdefault(go, set()).add(g)
    return go2genes


# -----------------------------
# Write GO term count file
# -----------------------------
def write_term_counts(outfile, go2genes, go_names=None):
    with open(outfile, "w") as out:
        for go in sorted(go2genes):
            genes = sorted(go2genes[go])
            name = go_names.get(go, go) if go_names else go
            out.write(f"{go}\t{name}\t{len(genes)}\t{','.join(genes)}\n")


# -----------------------------
# Precompute perm-index -> file list with a single directory scan
# -----------------------------
def build_perm_file_map(permdir, n_perms):
    perm_map = {i: [] for i in range(n_perms)}
    for fname in os.listdir(permdir):
        if "_perm_" not in fname:
            continue
        tail = fname.rsplit("_perm_", 1)[1]
        idx_str = tail.rsplit(".", 1)[0]
        try:
            i = int(idx_str)
        except ValueError:
            continue
        if i in perm_map:
            perm_map[i].append(os.path.join(permdir, fname))
    return perm_map


# -----------------------------
# MAIN
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--permdir", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--go", required=True)
    ap.add_argument("--n_perms", type=int, required=True)
    ap.add_argument("--term_outdir", required=True)
    ap.add_argument("--go_names", required=False, default=None,
                     help="Optional tab-separated GO_ID<TAB>Name file. If omitted, "
                          "the GO ID is written in both the ID and name columns "
                          "(same as the original script's behavior).")
    args = ap.parse_args()
    os.makedirs(args.term_outdir, exist_ok=True)

    print("loading genes")
    genes = load_genes(args.gff)
    print("genes:", len(genes))

    print("building gene index for fast overlap")
    gene_index = build_gene_index(genes)

    print("loading GO")
    gene2go = read_go(args.go)
    print("GO genes:", len(gene2go))

    go_names = read_go_names(args.go_names) if args.go_names else None

    tested_go_terms = build_tested_go_cache(genes, gene2go)
    print("tested GO terms:", len(tested_go_terms))

    # ---------------- REAL ----------------
    print("processing REAL")
    real_svs = read_intervals(args.real)
    real_hits = overlap(real_svs, gene_index)

    real_go_sets = fresh_tested_go_sets(tested_go_terms)
    real_hits_go_sets = go_gene_sets(real_hits, gene2go)
    for go in real_hits_go_sets:
        real_go_sets[go] = real_hits_go_sets[go]

    real_counts = {go: len(gset) for go, gset in real_go_sets.items()}

    write_term_counts(
        os.path.join(args.term_outdir, "../real_term_counts.tsv"),
        real_go_sets,
        go_names,
    )

    print("real SVs:", len(real_svs))
    print("real SVs with hits:", sum(len(h[4]) > 0 for h in real_hits))

    # ---- REAL gene tracking ----
    real_gene_set = set()
    for _, _, _, _, gene_set in real_hits:
        real_gene_set.update(gene_set)

    # ---------------- PERM ----------------
    print("processing PERM")

    perm_file_map = build_perm_file_map(args.permdir, args.n_perms)

    perm_dist = {}
    perm_gene_map = {}   # gene -> set of perm IDs where it appears

    for i in range(args.n_perms):

        svs = []
        for fpath in perm_file_map.get(i, []):
            svs.extend(read_intervals(fpath))

        hits = overlap(svs, gene_index)

        perm_go_sets = fresh_tested_go_sets(tested_go_terms)
        perm_hit_go_sets = go_gene_sets(hits, gene2go)
        for go in perm_hit_go_sets:
            perm_go_sets[go] = perm_hit_go_sets[go]

        counts = {go: len(gset) for go, gset in perm_go_sets.items()}

        write_term_counts(
            os.path.join(args.term_outdir, f"perm_{i}.tsv"),
            perm_go_sets,
            go_names,
        )

        # GO stats accumulation
        for go in real_counts:
            perm_dist.setdefault(go, []).append(counts.get(go, 0))

        # ---- PERM gene tracking ----
        perm_genes = set()
        for _, _, _, _, gene_set in hits:
            perm_genes.update(gene_set)

        for g in perm_genes:
            perm_gene_map.setdefault(g, set()).add(i)

        if (i + 1) % 10 == 0 or i == args.n_perms - 1:
            print(f"finished {i + 1}/{args.n_perms} perms")

    # ---------------- OUTPUT (GO) ----------------
    print("\nGO\tREAL\tMEAN_PERM\tENRICHMENT\tPVAL")

    for go, real_c in real_counts.items():

        perms = perm_dist.get(go, [0] * args.n_perms)
        mean_perm = np.mean(perms)

        if real_c == 0 and mean_perm == 0:
            enrichment = float("nan")   # no signal anywhere -- not "infinite"
        elif mean_perm > 0:
            enrichment = real_c / mean_perm
        else:
            enrichment = float("inf")   # real signal, zero permutation support

        pval = sum(p >= real_c for p in perms) / len(perms)

        print(go, real_c, mean_perm, enrichment, pval, sep="\t")

    # ---------------- OUTPUT (GENES) ----------------
    print("\n--- GENE OVERLAP SUMMARY ---")

    real_in_perms = {g for g in real_gene_set if g in perm_gene_map}

    print("REAL genes total:", len(real_gene_set))
    print("REAL genes also seen in at least one perm:", len(real_in_perms))

    print("\nGene\tPerm_indices")
    for g in sorted(real_in_perms):
        perms = sorted(perm_gene_map[g])
        print(g, ",".join(map(str, perms)), sep="\t")

    unique_real = real_gene_set - set(perm_gene_map.keys())

    print("\nREAL-only genes:", len(unique_real))
    for g in sorted(unique_real):
        print(g)


if __name__ == "__main__":
    main()
