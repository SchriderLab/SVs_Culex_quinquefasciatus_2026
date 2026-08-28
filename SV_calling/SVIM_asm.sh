#SVIM-asm calling commands run for each sample
minimap2 -ax asm5 --cs -r2k -t 8 culex.qui.gcf.fasta /proj/dschridelab/remi/ANNA_RELAB/01_2026_FINAL2/01_GENOME_ASS/LZV1_GENOME_FINAL2.fasta > LZV1.sam
samtools sort -m4G -@4 -o LZV1.sorted.bam LZV1.sam
samtools index LZV1.sorted.bam
svim-asm haploid ./LZV1 LZV1.sorted.bam culex.qui.gcf.fasta
