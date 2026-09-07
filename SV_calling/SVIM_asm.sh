#SVIM-asm calling commands run for each sample
minimap2 -ax asm5 --cs -r2k -t 8 GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna [sample]_GENOME_FINAL2.fasta > [sample].sam
samtools sort -m4G -@4 -o [sample].sorted.bam [sample].sam
samtools index [sample].sorted.bam
svim-asm haploid ./[sample] [sample].sorted.bam GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna
