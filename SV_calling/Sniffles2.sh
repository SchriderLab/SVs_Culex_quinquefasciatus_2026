#Sniffles calling commands run for each sample
sniffles -i culqui.gcf.LZV1.NEW.bam -v LZV1_sniffles.vcf  --max-splits-kb 0.1 --max-splits-base 2 --long-del-length 10000 --minsupport 8 \
--reference GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna
