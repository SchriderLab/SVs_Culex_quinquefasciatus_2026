#pbsv calling commands run for each sample
pbsv discover culqui.gcf.LZV1.NEW.bam LZV1.svsig.gz --hifi
pbsv call --log-level INFO --hifi -j 20 GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fasta LZV1.svsig.gz LZV1_pbsv.vcf
