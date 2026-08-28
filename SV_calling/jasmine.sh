# To comply with the formatting requirements for Jasmine v1.1.5, we ran the following steps:
# Added arbitrary strand values & added strands to vcf header (run jasmine with --ignore_strand option) (02_INT_FILES/02_JASMINE_FORMATTING/*_strands.vcf)
# Removed all instances of "." from the VCF to avoid java.lang.NumberFormatException bug when running Jasmine (02_INT_FILES/02_JASMINE_FORMATTING/*_reformatted.vcf)
# Removed the ".1" from chrom names and variant IDs and changed all qual scores to 100 (arbitrary). Qual scores are all set to "." after running combiSV, and jasmine doesn't take qual scores into account (02_INT_FILES/02_JASMINE_FORMATTING/*_qual_reformatted.vcf)
# Ensured that final preprocessed vcf files were tab-separated instead of space-separated before running Jasmine (02_INT_FILES/02_JASMINE_FORMATTING/*_qual_reformatted_tab.vcf)
# Added ".1" back to chromosomes in Jasmine output file before splitting samples and running Sniffles (02_INT_FILES/03_JASMINE_RUN/jasmine_filt_pregt_chromfixed.vcf)

#1 - reformat VCFs
# 1.1 - add strands & matching VCF header
awk 'BEGIN{FS=OFS="\t"}
/^##/ {print; next}
/^#CHROM/ {
    print "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strands\">"
    print
    next
}
{
    $8 = $8 ";STRANDS=+-"
    print
}' LZV1_combisv_2callers_sorted.vcf > LZV1_strands.vcf

# 1.2 - remove ".1" from chromosome names and "." from SV ids
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    sub(/\..*$/, "", $1)
    sub(/^id\./, "", $3)
    print
}' LZV1_strands.vcf > LZV1_reformatted.vcf

# 1.3 - change qual scores from "." to 100
grep "#" LZV1_reformatted.vcf > LZV1_qual_reformatted.vcf
grep -v "#" LZV1_reformatted.vcf | awk -i inplace '{$6=100; print}' >> LZV1_qual_reformatted.vcf

# 1.4 - for GTs pulled from SVIM - get rid of "." in GT qualities
sed -i 's/\./0/g' LZV1_qual_reformatted.vcf

#2 - run Jasmine
# 2.1 - combine combisv 2+ caller VCFs with Jasmine
jasmine file_list=combisv_jasmine_samples.txt out_file=combisv_jasmine_merged.vcf --ignore_strand --output_genotypes
#file_list = paths to the combiSV VCFs (filtered for 2+ caller support and sorted)

# add ".1" back to chromosome names after running Jasmine & before splitting samples (needed for Sniffles genotyping/downstream analysis)
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    $1 = $1 ".1"
    print
}' combisv_jasmine_merged.vcf > jasmine_filt_pregt_chromfixed.vcf


# 2.2 - split Jasmine output into single-sample vcfs to run Sniffles genotyping
bcftools view -s LZV1 jasmine_filt_pregt_chromfixed.vcf > LZV1_jasmine.vcf

# 2.3 - Sniffles genotyping/force-calling step - run on each sample
sniffles --input culqui.gcf.LZV1.NEW.bam --genotype-vcf LZV1_jasmine.vcf --vcf LZV1_gt.vcf

# 2.4 - combine genotyped samples (must gzip & index genotyped VCFs first)
bcftools view LZV1_gt.vcf -Oz -o LZV1_gt.vcf.gz
bcftools index LZV1_gt.vcf.gz
bcftools merge *_gt.vcf.gz > jasmine_merged.vcf

# 2.5 - remove any variants where all samples were genotyped as homozygous reference (0/0)
grep "#" jasmine_merged.vcf > jasmine_filt_no0.vcf
grep -v "AC=0" jasmine_merged.vcf >> jasmine_filt_no0.vcf

# 2.6 - generate vcf with no missing data for PCA and SFS
bcftools view -e 'GT[*]="mis"' jasmine_filt_no0.vcf -o jasmine_filt_nomissing_no0.vcf
