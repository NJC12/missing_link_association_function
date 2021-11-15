for chr in {1..22}; do
cat > ~/tmp/calc.ukbb.mafs.chr${chr}.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/calc.ukbb.mafs.chr${chr}.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/calc.ukbb.mafs.chr${chr}.sh.o%j
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=50G
#SBATCH -p medium
#SBATCH -t 5-00:00:00

mkdir -p /n/scratch3/users/n/[username]/gwas/results/maf

module load plink2/2.0

# Need to include --covar-variance-standardize or else I get an error.
# The intermediate files that plink produces are large, so I need to write to scratch to keep my jobs from dying for lack of room in my home directory.

plink2 \\
--bgen      /n/groups/price/UKBiobank/bgen_MAF001_500K_v3/UKB_MAF0.001_v3.${chr}.bgen \\
--sample    /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample \\
--remove    /n/groups/price/UKBiobank/sampleQC/remove.nonStringentBritish.FID_IID.txt \\
        /n/groups/price/UKBiobank/sampleQC/remove.related.FID_IID.txt \\
        /n/groups/price/UKBiobank/download_500K/w14048_20200204.FID_IID.txt \\
--exclude   /n/groups/price/UKBiobank/snpQC/autosome_maf_lt_0.001.txt \\
        /n/groups/price/UKBiobank/snpQC/autosome_missing_gt_0.1.txt \\
        /home/[username]/pleiotropy_o2/gwas/data/snps/cholesterol/coding.snps.for.cholesterol.tsv \\
--threads   1 \\
--freq      \\
--out       /n/scratch3/users/n/[username]/gwas/results/maf/ukbb.maf.chr${chr}

EOF
done
