#1. we followed the 0. Quality Control and Read Mapping and 1. Variant Calling and Filtering workflows from Haba et al. (2025), found in https://github.com/YukiHaba/PipPop_molestus_origin/tree/main/scripts
# with the following modifications:
  # --include 'F_MISSING < 0.25 & MQ > 40 & QUAL > 30' \ in 1. variant_calling.sh was changed to --include 'F_MISSING < 0.5 & MQ > 40 & QUAL > 30' \
  # we ran the scripts only on the samples identified as C. quinquefasciatus in the original study

#2. subset VCFs by population (South Africa, Madagascar, Cameroon/Gabon)
bcftools view -S safrica.txt good_biallelic_snps_05.rm.combined.accessible.vcf.gz -Oz -o snps_final_southafrica.vcf.gz

#3. run vcf2sf.py (https://github.com/SchriderLab/timesweeper-experiments/blob/main/scripts/comp_methods/vcf2sf.py) - generates the "FreqFile" for each chromosome to use for SweepFinder

#4. SweepFinder commands:
# "CombinedFreqFile" = the 3 "FreqFiles" for each chromosome in a population concatenated together

# run the following for each population to generate SFS:
SweepFinder2 -f CombinedFreqFile SpectFile

# run the following for each population. "OutFile" = CLR results
SweepFinder2 -lg 1000 FreqFile SpectFile OutFile

