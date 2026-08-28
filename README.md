# SVs_Culex_quinquefasciatus_2026

## Genome assembly and annotation

## Species confirmation
`species_ID/mitochondrial_tree.sh`: commands for generating the mitochondrial tree in Supplemental Figure 1\

`species_ID/nuclear_tree.sh`: commands for generating the nuclear genome tree in Figure 1


## Structural variant calling
`SV_calling` contains the commands used to call SVs in each sample for the four individual callers used here, as well as:\

`combiSV.sh`: combiSV commands used to merge evidence from the four callers for each sample\

`combiSV2.3_modified.pl`: our modified combiSV script where the logic requiring minimum read support was removed for the SVIM vcf\

`jasmine.sh`: commands for combining combiSV results into a multi-sample VCF using Jasmine, individually genotyping these merged variants in each sample using Sniffles2, and generating the final multi-sample full set VCF\

`generate_final_files.sh`: commands used to generate our six additional VCFs: the high-impact VCF, outlier VCF, high-impact outlier VCF, LZV-enriched VCF, PR-enriched VCF, and LZV/PR enriched VCF (and commands to filter each VCF for SVs overlapping IR genes)\

`vcf2bed.sh`: script used in `generate_final_files.sh` to extract SV coordinates


## SNP calling
`SNP_calling/call_snps.sh`: commands for joint calling SNPs in our *C. quinquefasciatus* samples and annotating with SnpEff


## SweepFinder analysis
`sweepfinder/sweepfinder.sh`: commands for running SweepFinder analysis described in Supplemental XXX. Detailed methods for QC, read mapping, and variant calling can be found [here](https://github.com/YukiHaba/PipPop_molestus_origin/tree/main) 
