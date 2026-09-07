#SVIM-asm calling commands run for each sample
minimap2 -ax asm5 --cs -r2k -t 8 GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna LZV1_GENOME_FINAL2.fasta > LZV1.sam
samtools sort -m4G -@4 -o LZV1.sorted.bam LZV1.sam
samtools index LZV1.sorted.bam
svim-asm haploid ./LZV1 LZV1.sorted.bam GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna
