# Structural Variants in *Culex quinquefasciatus*

## Genome assembly 
`Genome_Assembly/assembly.sh`: Assemble all of the *C. quinquefasciatus* genomes


## Species confirmation
`species_ID/mitochondrial_tree.sh`: Generate the mitochondrial tree

`species_ID/nuclear_tree.sh`: Generate the nuclear genome tree


## Structural variant calling
`SV_calling/`: Call SVs in each sample for the four individual callers used here

`SV_calling/combiSV.sh`: Merge evidence from the four callers for each sample  

`SV_calling/combiSV2.3_modified.pl`: Modified combiSV script where the logic requiring minimum read support was removed for the SVIM vcf  

`SV_calling/jasmine.sh`: Combine combiSV results into a multi-sample VCF using Jasmine, individually genotyping these merged variants in each sample using Sniffles2, and generating the final multi-sample full set VCF  

`SV_calling/generate_final_files.sh`: Generate our six additional VCFs: the high-impact VCF, outlier VCF, high-impact outlier VCF, LZV-enriched VCF, PR-enriched VCF, and LZV/PR enriched VCF (and commands to filter each VCF for SVs overlapping IR genes)  

`SV_calling/vcf2bed.sh`: script used in `generate_final_files.sh` to extract SV coordinates


## SNP calling
`SNP_calling/call_snps.sh`: commands for joint calling SNPs in our *C. quinquefasciatus* samples and annotating with SnpEff


## SweepFinder analysis
`sweepfinder/sweepfinder.sh`: Run SweepFinder as described in the Supplemental Material. Detailed methods for QC, read mapping, and variant calling can be found [here](https://github.com/YukiHaba/PipPop_molestus_origin/tree/main) 

## Enrichment analysis
`Enrichment/01_EXTRACT_MASTER_SV.sh`
`Enrichment/02_PERMUTE_MASTER_SV.sh`
`Enrichment/03_FILTER_PERMS_BY_DATASET.py`
`Enrichment/04_FILTER_REAL_BYDATSET.py`
`Enrichment/05_GET_TERM_COUNTS.py`
`Enrichment/06_GET_PERMUTAT_QVALS_GENCOUNTS.py`
`Enrichment/07_APPEND_GO_NAMES.py`
`Enrichment/09_DO_PERM_ON_GENELIST.py`

