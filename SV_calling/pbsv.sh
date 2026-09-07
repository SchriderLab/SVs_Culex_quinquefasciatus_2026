#pbsv calling commands run for each sample
#generate bams first
pbmm2 align GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fasta LZV1.fastq culqui.gcf.LZV1.NEW.bam --sort --preset CCS --sample LZV1 --rg '@RG\tID:LZV1'


#run pbsv
pbsv discover culqui.gcf.LZV1.NEW.bam LZV1.svsig.gz --hifi
pbsv call --log-level INFO --hifi -j 20 GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fasta LZV1.svsig.gz LZV1_pbsv.vcf
