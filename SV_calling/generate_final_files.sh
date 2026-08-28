#files we generated using the full set VCF:

#High-impact VCF
#1 - run SnpEff on full set VCF

# 1.1 - build genome file. need the following files from NCBI:
# /snpEff/./data/GCF_015732765.1/cds.fa
# /snpEff/./data/GCF_015732765.1/protein.fa
# /snpEff/./data/GCF_015732765.1/genes.gtf
# /snpEff/./data/genomes/GCF_015732765.1.fa
# genome must be added to config file
java -jar snpEff.jar build -gtf22 -v GCF_015732765.1

# 1.2 - run SnpEff
java -jar snpEff.jar ann GCF_015732765.1 jasmine_filt_no0_under100kb.vcf \
-csvStats jasmine_filt_stats.csv -v > jasmine_filt_snpeff.vcf

# 1.3 - filter for high-impact SVs
grep "#" jasmine_filt_snpeff.vcf > highimpact.vcf
grep "HIGH" jasmine_filt_snpeff.vcf >> highimpact.vcf


#Outlier VCF
#1. get positions above 99th percentile for CLR value from sweepfinder output for each chromosome
sf_data <- read.table("/Users/annatrotter/Desktop/lab/sweepfinder/sweepfinder_cg_chr3_merged.txt", header=TRUE)
print(quantile(sf_data$LR,0.99))
awk -F '\t' 'NR==1 || $2 > “X”’ sweepfinder_cg_chr3_merged.txt > cg_chr3_99.txt

#2. run sweepfinder_ranges.py for each pop/chromosome to generate peak region coordinates. add chrom names and concatenate/sort files to create sweepfinder_ranges_all.bed

#3. convert SV coordinates from the full set VCF to BED format using vcf2bed.sh

#4. find intersections between sweepfinder ranges and SVs
bedtools intersect -wa -wb -a jasmine_filt_no0_under100kb.bed -b sweepfinder_ranges_all.bed > outlier_SVs.txt

#5. extract SV IDs from outlier_SVs.txt and filter the main VCF
(
          	grep '^#' jasmine_filt_no0_under100kb.vcf
                grep -Fwf 01_INT_FILES/outlier_ids.txt jasmine_filt_no0_under100kb.vcf
        ) > outlier.vcf


#High-impact outlier VCF
#1. filter outlier VCF for "HIGH"
grep "#" outlier.vcf > highimpact_outlier.vcf
grep "HIGH" outlier.vcf >> highimpact_outlier.vcf


#LZV-enriched VCF

# TODO: REMI

#PR-enriched VCF

# TODO: REMI

#LZV/PR-enriched VCF

# TODO: REMI

#we filtered each VCF for IR gene overlap using the following commands:
#note: the IR_genes.bed file created from A. gambiae orthologs info - chrom, start, end (see Supplemental Data)

#1 - get SVs overlapping IR genes
bedtools intersect -wa -wb -a jasmine_filt_no0_under100kb.bed -b IR_genes.bed > jasmine_filt_overlap_IRgenes.bed
#IR gene IDs were extracted from this BED file & duplicates removed > jasmine_overlap_IRgene_ids.txt

#2- filter for SVs overlapping IR genes
(
                grep '^#' jasmine_filt_no0_under100kb.vcf
                grep -Fwf jasmine_overlap_IRgene_ids.txt jasmine_filt_no0_under100kb.vcf
        ) > jasmine_overlapping_IRgenes.vcf

#3 - for the high-impact VCF, we filtered using the following command to isolate variants exerting high-impact effects specifically on IR genes
(
                grep '^#' 01_INT_FILES/jasmine_overlapping_IRgenes_snpeff.vcf
                grep -Ff <(sed 's/^/HIGH|LOC/' 01_INT_FILES/IR_ids.txt) 01_INT_FILES/jasmine_overlapping_IRgenes_snpeff.vcf
        ) > jasmine_overlapping_IRgenes_high_snpeff.vcf
