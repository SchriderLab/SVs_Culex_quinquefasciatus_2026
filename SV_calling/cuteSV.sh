#cuteSV calling command run for each sample - default HiFi parameters + --max-split-parts filter

cuteSV ./culqui.gcf.LZV1.NEW.bam ./GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna LZV1_cutesv.vcf ./LZV1 \
--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max-cluster_bias_DEL 1000 \
 --diff_ratio_merging_DEL 0.5 -s 5 --max_split_parts 4 --genotype
