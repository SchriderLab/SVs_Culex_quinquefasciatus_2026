#1 - run mitogenome assembly for each sample using MitoHiFi
singularity exec --bind /work/users/a/n/annawt/SPECIESID_NEW/mitohifi/ --bind /proj/dschridelab/remi/ANNA_RELAB:/proj/dschridelab/remi/ANNA-RELAB docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py \
-r /proj/dschridelab/remi/ANNA_RELAB/01_2026_FINAL2_RNK/00_FASTQ/LZV1.fastq.gz -f NC_014574.1.fasta -g NC_014574.1.gb -t 72

#2 - rotate mitogenomes to a common starting point using Rotate
./rotate/rotate -s TCGCGACAATG -m 1 culex_aedes_mito.fa > culex_aedes_mito.rotated.fa

#3 - perform MAFFT alignment (done in UGENE)

#4 - run IQ-TREE on this alignment
iqtree2 -s culex_aedes_mito_alignment.fa -B 1000
