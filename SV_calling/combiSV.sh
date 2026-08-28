#combiSV command to merge evidence from the 4 callers run on each sample
#the combiSV script was modified to remove read support logic for SVIM, since we used assembly-based SVIM-asm

#1 - run combiSV for each sample
#perl combiSV2.3_modified.pl -cutesv LZV1_cutesv.vcf -sniffles LZV1_sniffles.vcf -svim LZV1_svim.vcf -pbsv LZV1_pbsv.vcf -o LZV1_combisv.vcf

#2 - filter for SVs supported by 2+ callers
grep "#" LZV1_combisv.vcf > LZV1_combisv_2callers.vcf
grep 'SVCALLERS=Sniffles,pbsv,cutesv,SVIM' >>  LZV1_combisv_2callers.vcf
grep -E 'SVCALLERS=(Sniffles,pbsv,cutesv|Sniffles,pbsv,SVIM|Sniffles,cutesv,SVIM|pbsv,cutesv,SVIM)(     |$)' >> LZV1_combisv_2callers.vcf
grep -E 'SVCALLERS=(Sniffles,pbsv|Sniffles,cutesv|Sniffles,SVIM|pbsv,cutesv|pbsv,SVIM|cutesv,SVIM)(     |$)' >> LZV1_combisv_2callers.vcf

#3 - sort the VCF
bcftools sort LZV1_combisv_2callers.vcf -o LZV1_combisv_2callers_sorted.vcf
