import numpy as np

vcf_file = "jasmine_filt_no0_under100kb.vcf"

pop1 = ["LZV1", "LZV3", "LZV4", "LZV5"]
pop2 = ["PR1", "PR4", "PR9"]

sv_records = []
header_lines = []


def allele_freq(genotypes, indices):

    alt = 0
    total = 0

    for i in indices:

        gt = genotypes[i]

        if "." in gt:
            continue

        alleles = gt.replace("|", "/").split("/")

        alt += sum(int(a) for a in alleles)
        total += 2

    return alt / total if total else np.nan


#read VCF and calc AF differences
with open(vcf_file) as f:

    for line in f:

        if line.startswith("#"):

            header_lines.append(line)

            if line.startswith("#CHROM"):

                header = line.strip().split("\t")
                samples = header[9:]

                pop1_idx = [samples.index(x) for x in pop1]
                pop2_idx = [samples.index(x) for x in pop2]


        else:

            fields = line.strip().split("\t")

            genotypes = []

            for sample in fields[9:]:
                gt = sample.split(":")[0]
                genotypes.append(gt)


            af_lzv = allele_freq(genotypes, pop1_idx)
            af_pr = allele_freq(genotypes, pop2_idx)


            if not np.isnan(af_lzv) and not np.isnan(af_pr):

                # Positive = LZV enriched
                # Negative = PR enriched
                diff = af_lzv - af_pr

                sv_records.append((diff, line))


# Calculate percentile cutoffs
differences = np.array([x[0] for x in sv_records])

bottom_cutoff = np.percentile(differences, 5)
top_cutoff = np.percentile(differences, 95)


print("Total SVs analysed:", len(sv_records))
print("PR-enriched cutoff (5th percentile):", bottom_cutoff)
print("LZV-enriched cutoff (95th percentile):", top_cutoff)


# Output filenames
combined_output = "SVs_extreme_AF_difference_top_bottom5.vcf"
pr_output = "PR_enriched_SVs_AF_difference_bottom5.vcf"
lzv_output = "LZV_enriched_SVs_AF_difference_top5.vcf"


# Write combined extreme 10%
with open(combined_output, "w") as out:

    for h in header_lines:
        out.write(h)

    for diff, line in sv_records:

        if diff <= bottom_cutoff or diff >= top_cutoff:
            out.write(line)


# Write PR enriched bottom 5%
with open(pr_output, "w") as out:

    for h in header_lines:
        out.write(h)

    for diff, line in sv_records:

        if diff <= bottom_cutoff:
            out.write(line)


# Write LZV enriched top 5%
with open(lzv_output, "w") as out:

    for h in header_lines:
        out.write(h)

    for diff, line in sv_records:

        if diff >= top_cutoff:
            out.write(line)


print("\nOutput files:")
print(combined_output)
print(pr_output)
print(lzv_output)
