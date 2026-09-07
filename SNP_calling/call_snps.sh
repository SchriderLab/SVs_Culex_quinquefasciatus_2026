#joint call SNPs using bcftools - broken into regions to parallelize (this is an example for one region)
bcftools mpileup --config pacbio-ccs -Ou -f GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \
culqui.gcf.LZV1.NEW.bam culqui.gcf.LZV3.NEW.bam culqui.gcf.LZV4.NEW.bam \
culqui.gcf.LZV5.NEW.bam culqui.gcf.PR1.NEW.bam culqui.gcf.PR4.NEW.bam culqui.gcf.PR9.NEW.bam \
--max-depth 1000 -r NC_051863.1:1-50000000 | bcftools call -mv -Ob -o all_calls_chr3:start-50.bcf


#concatenate regions into one VCF and sort
bcftools concat all_calls_chr1:start-50.bcf \
all_calls_chr1:50-100.bcf all_calls_chr1:100-end.bcf \
all_calls_chr2:start-50.bcf all_calls_chr2:50-100.bcf all_calls_chr2:100-150.bcf \
all_calls_chr2:150-200.bcf all_calls_chr2:200-end.bcf \
all_calls_chr3:start-50.bcf all_calls_chr3:50-100.bcf all_calls_chr3:100-150.bcf all_calls_chr3:150-end.bcf  -O z -o joint_snps_all.vcf.gz

bcftools sort joint_snps_all.vcf.gz -O z -o joint_snps_all.sorted.vcf.gz


#filter for QUAL >= 20
bcftools view -i 'QUAL>=20' joint_snps_all.sorted.vcf.gz -O z -o joint_snps_all.sorted.Q20.vcf.gz


#run SnpEff (see SV_calling/generate_final_files.sh for more detailed description of setup)
java -jar snpEff/snpEff.jar ann GCF_015732765.1 joint_snps_all.sorted.Q20.vcf.gz -csvStats joint_snps_all.sorted.Q20.csv -v \
| bcftools view -O z -o joint_snps_all.sorted.Q20.snpeff.vcf.gz
