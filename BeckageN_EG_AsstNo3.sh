#!/bin/bash
CHR="Chr15"

echo $CHR  # Does it look right?

cd /data/project_data/PopGenomics/shared/${CHR}_sweeps/LAI/

# Get phenotype data from meta file for the admixed samples:

tail -n +2 /data/project_data/PopGenomics/climDat.txt | \
grep -w -f Admixed.Inds - >Admixed.pheno

printf "#FID\tIID\tmean_finalFreeze\tmean_cGDDfreeze\tmed_DD0\n%.0s" >Plink/phenoHW.forPlink

cat Admixed.pheno >>Plink/phenoHW.forPlink

# Get K=2 ADMIX to use as covariate; Need to use tail -n +2 to skip the header in the metadata file before pasting to the Q file

tail -n +2 /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1 | \
paste - /data/project_data/PopGenomics/shared/poplar_hybrids.LDpruned.5.Q | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f 1,5 >Plink/Admixed_KBals

cat /data/project_data/PopGenomics/Combined_Transect_Sampling_Data_2020.txt | \
cut -f1-2 | \
grep -w -f Admixed.Inds - | \
sed 's/\s/\t/g' | \
cut -f2 >Plink/Transect

# Create the cov file with KBals

printf "#FID\tIID\tK2Q\n%.0s" >Plink/covHW.forPlink
paste Plink/Transect Plink/Admixed_KBals >>Plink/covHW.forPlink

cd Plink/

plink2 --bfile ${CHR}_Admixed_FAI \
--pheno phenoHW.forPlink \
--covar covHW.forPlink \
--glm omit-ref