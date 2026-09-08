# Structural Variants in *Culex quinquefasciatus*

## Genome assembly 
`Genome_Assembly/assembly.sh`: Assemble all of the *C. quinquefasciatus* genomes


## Species confirmation
`species_ID/mitochondrial_tree.sh`: Generate the mitochondrial tree

`species_ID/nuclear_tree.sh`: Generate the nuclear genome tree


## Structural variant calling
`SV_calling/Sniffles2.sh`, `SV_calling/cuteSV.sh`, `SV_calling/pbsv.sh`, SV_calling/SVIM_asm.sh`: Call SVs in each sample with the four individual callers used here

`SV_calling/combiSV.sh`: Merge evidence from the four callers for each sample  

`SV_calling/combiSV2.3_modified.pl`: Modified combiSV script where the logic requiring minimum read support was removed for the SVIM vcf  

`SV_calling/jasmine.sh`: Combine combiSV results into a multi-sample VCF using Jasmine, individually genotyping these merged variants in each sample using Sniffles2, and generating the final multi-sample full set VCF  

`SV_calling/generate_final_files.sh`: Generate our six additional VCFs: the high-impact VCF, outlier VCF, high-impact outlier VCF, LZV-enriched VCF, PR-enriched VCF, and LZV/PR enriched VCF (and commands to filter each VCF for SVs overlapping IR genes)  

`SV_calling/vcf2bed.sh`: Extract SV coordinates from VCF as BED file (used in `generate_final_files.sh`) 

`SV_calling/sweepfinder_ranges.py`: Extract ranges of CLR peaks based on SweepFinder analysis (used to generate our outlier VCFs)  

`SV_calling/extract_AF_VCFs.py`: Generate LZV, PR, and LZV/PR-enriched VCFs


## SNP calling
`SNP_calling/call_snps.sh`: Commands for joint calling SNPs in our *C. quinquefasciatus* samples and annotating with SnpEff


## SweepFinder analysis
`sweepfinder/sweepfinder.sh`: Run SweepFinder as described in the Supplemental Material. Detailed methods for QC, read mapping, and variant calling can be found [here](https://github.com/YukiHaba/PipPop_molestus_origin/tree/main) 

## Enrichment analysis
`Enrichment/01_EXTRACT_MASTER_SV.sh`: Extract all the SV info (SV ID, chrom, position start, position finish) from jasmine_filt_no0_under100kb.vcf (full SV list). 

`Enrichment/02_PERMUTE_MASTER_SV.sh`: Permute the master SV file. This circular shift uses a single random offset per chrom per permutation.  

`Enrichment/03_FILTER_PERMS_BY_DATASET.py`: Filter the master permuted SV files by SVs present in VCF of interest.  

`Enrichment/04_FILTER_REAL_BYDATSET.py`: Filter the master SV interval file (not permuted) by VCF of interest.  

`Enrichment/05_GET_TERM_COUNTS.py`: Test whether genes affected by SVs are unusually enriched for particular Gene Ontology (GO) functions by comparing the real data against permutations  

`Enrichment/06_GET_PERMUTAT_QVALS_GENCOUNTS.py`: Calc p-values, enrichment scores, and q-values (FDR-adjusted significance) for each GO term by comparing the real gene counts with permutation gene counts. 

`Enrichment/07_APPEND_GO_NAMES.py`: Takes GO annotation file and the enrichment results file and writes a tab-sep file where each GO ID from the enrichment results is annotated with all relevant info.  

`Enrichment/08_DO_PERM_ON_GENELIST.py`: Test whether IR genes overlap SV regions more often than expected by chance.

