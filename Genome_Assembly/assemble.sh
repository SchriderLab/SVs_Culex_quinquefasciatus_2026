## Code below was used to assemble the LZV3 genome, note that some of the exact parameters varied based on which individual we were assembling.


module load samtools
module load minimap2/2.26

conda activate /work/users/r/k/rketchum/software/purge_haplotigs/
conda activate /work/users/r/k/rketchum/software/kat/
conda activate /work/users/r/k/rketchum/software/ragtag/
conda activate /work/users/r/k/rketchum/software/purge_dups/



# Assemble with hifiasm then convert gfa to fa

./hifiasm/hifiasm -o [output] -t 40 [sample].fastq.gz
awk '/^S/{print ">"$2;print $3}' [sample].bp.p_ctg.gfa > [sample].bp.p_ctg.fa


# Check the kat kmer plots to take a look at how well the genomes have been purged

kat comp -t 16 -o [sample].KAT [sample].fastq.gz [sample].bp.p_ctg.fa

# To purge

minimap2 -t16 -ax map-pb [sample].bp.p_ctg.fa [sample].fastq.gz --secondary=no | samtools sort -m 250G -o [sample].p_ctg.bam -T tmp.ali11

samtools index [sample].p_ctg.bam

purge_haplotigs hist -b [sample].p_ctg.bam  -g [sample].bp.p_ctg.fa -t 16

purge_haplotigs contigcov -i [sample].p_ctg.bam.200.gencov -l 15 -m 90  -h 170  -o [sample].cov

purge_haplotigs purge -t 16 -g [sample].bp.p_ctg.fa -c [sample].cov -a 50 -o [sample].purge

purge_haplotigs clip -t 16 -l 5000 -p [sample].purge.fasta -h [sample].purge.haplotigs.fasta -o [sample].purge.clip


# Check kat kmer plots for newly purged genome assembly

kat comp -t 16 -o [sample].KAT.CLIP [sample].fastq.gz [sample].purge.clip.fasta

# Purge again

run_purge_dups.py config.json ./purge_dups/bin/ [sample]_PUR_CLIP -p bash

# Check kmer again 

kat comp -t 16 -o [sample].PURDUPS [sample].fastq [sample].purge.clip.purged.fa


# Rag Tag to stitch it all together

ragtag.py correct \
./GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \
-R [sample].fastq -T corr [sample].purge.clip.purged.fa

ragtag.py scaffold \
./GCF_015732765.1_VPISU_Cqui_1.0_pri_paternal_genomic.fna \
ragtag_output2/ragtag2.correct.fasta

# remove short sequences

seqkit seq -m 1000 ragtag_output2/ragtag2.scaffold.fasta > [sample]_FINAL2.fasta
