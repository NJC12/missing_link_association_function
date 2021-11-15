#!/bin/bash

### Create phenotype and covariate files

cat > ~/tmp/create.height.pheno.file.sh <<'EOF'
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/create.height.pheno.file.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/create.height.pheno.file.sh.o%j
#SBATCH --mem=15G
#SBATCH -p short
#SBATCH -t 5



# MAKE SURE THIS `NR > 1` IS CORRECT IN JUST SKIPPING THE PHENO NAME.

awk '
NR == FNR && NR > 1 {
key[$1]=$2
}
NR != FNR && FNR == 1 {
print "FID\tIID\t" $10
}
NR != FNR && FNR > 1 {
print key[$1] "\t" key[$1] "\t" $10
}' \
/n/groups/price/UKBiobank/app10438assoc/ukb4989.IDmap.tab \
/n/groups/price/UKBiobank/app10438assoc/ukb4986.tab \
> /home/[username]/pleiotropy_o2/gwas/data/height.14048_ids.tab

awk '
{print $1 "\t" $2 "\t" $20 "\t" $21 "\t" $36 "\t" $73 "\t" $74 "\t" $75 "\t" $76 "\t" $77 "\t" $78 "\t" $79 "\t" $80 "\t" $81 "\t" $82}
' \
<(zcat /n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz) \
> /home/[username]/pleiotropy_o2/gwas/data/height.covariates.14048_ids.tab
EOF

scp ~/tmp/create.height.pheno.file.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/create.height.pheno.file.sh
'

### Extract height coding SNPs

mem_array=( 1000 256 256 256 256 256 256 256 256 150 150 150 150 150 150 150 150 100 150 100 100 100 100 )
time_array=( 1000 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4 2 6 2 2 2 2 )

for chr in {1..1}; do

cat > ~/tmp/create.height.pheno.file.chr$chr.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/create.height.pheno.file.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/create.height.pheno.file.sh.o%j
#SBATCH --mem=${mem_array[$chr]}G
#SBATCH -p short
#SBATCH -t ${time_array[$chr]}:00:00



module load plink2/2.0

plink2 \\
--bgen       /n/groups/price/UKBiobank/bgen_MAF001_500K_v3/UKB_MAF0.001_v3.$chr.bgen \\
--sample     /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample \\
--remove     /n/groups/price/UKBiobank/sampleQC/remove.nonStringentBritish.FID_IID.txt \\
           /n/groups/price/UKBiobank/sampleQC/remove.related.FID_IID.txt \\
           /n/groups/price/UKBiobank/download_500K/w14048_20200204.FID_IID.txt \\
--extract    /home/[username]/pleiotropy_o2/gwas/data/snps/height/coding.snps.for.height.tsv \\
--export     vcf \\
--out        /n/scratch3/users/n/[username]/gwas/data/snps/height/ukbb.height.coding.snps.chr$chr

EOF

done

scp ~/tmp/create.height.pheno.file.chr*.sh [username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..1}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/create.height.pheno.file.chr$chr.sh
done
'

### Transpose VCFs

for chr in {1..1}; do
cat > ~/tmp/run.transpose.vcfs.chr$chr.py <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.transpose.height.vcfs.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.transpose.height.vcfs.sh.o%j
#SBATCH --mem=25G
#SBATCH -p short
#SBATCH -t 00:30:00

module load gcc/6.2.0
module load python/3.7.4

python3.7 /home/[username]/pleiotropy_o2/gwas/tmp_src/transpose.vcfs.py \
/n/scratch3/users/n/[username]/gwas/data/snps/height/ukbb.height.coding.snps.chr$chr.vcf \
/n/scratch3/users/n/[username]/gwas/data/snps/height/transposed.ukbb.height.coding.snps.chr$chr.vcf

EOF
done

scp \
~/tmp/transpose.vcfs.py \
~/tmp/run.transpose.vcfs.chr*.py \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..1}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.transpose.vcfs.chr$chr.py
done
'

### Combine regular covariates with coding SNPs

for chr in {1..1}; do
cat > ~/tmp/run.combine.height.covariates.chr$chr.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.combine.height.covariates.chr$chr.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.combine.height.covariates.chr$chr.sh.o%j
#SBATCH --mem=5G
#SBATCH -p short
#SBATCH -t 00:10:00

mkdir -p /n/scratch3/users/n/[username]/gwas/data/snps/height
cd /home/[username]/pleiotropy_o2/gwas/data

awk '
NR == FNR {
covs[\$1]=\$0
}
NR != FNR {
if (\$1 in covs) {
  str=covs[\$1];
  for (i=3; i<=NF; i++) {
    str=str"\t"\$i
  }
}
print str
}' \\
/home/[username]/pleiotropy_o2/gwas/data/height.covariates.14048_ids.tab \\
/n/scratch3/users/n/[username]/gwas/data/snps/height/transposed.ukbb.height.coding.snps.chr$chr.vcf | \\
sort | uniq \\
> /n/scratch3/users/n/[username]/gwas/data/snps/height/height.covariates.w.snps.chr$chr.tsv

EOF
done

scp \
~/tmp/run.combine.height.covariates.chr*.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..1}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.combine.height.covariates.chr$chr.sh
done
'

### Run the GWAS

for chr in {1..22}; do
cat > ~/tmp/run.height.gwas.chr${chr}.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.height.gwas.chr${chr}.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.height.gwas.chr${chr}.sh.o%j
#SBATCH --mem=200G
#SBATCH -p long
#SBATCH -t 10-00:00:00

mkdir -p /n/scratch3/users/n/[username]/gwas/results/height

module load plink2/2.0

# Need to include --covar-variance-standardize or else I get an error.
# The intermediate files that plink produces are large, so I need to write to scratch to keep my jobs from dying for lack of room in my home directory.

plink2 \\
--bgen      /n/groups/price/UKBiobank/download_500K/ukb_imp_chr${chr}_v3.bgen \\
--sample    /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample \\
--remove    /n/groups/price/UKBiobank/sampleQC/remove.nonStringentBritish.FID_IID.txt \\
          /n/groups/price/UKBiobank/sampleQC/remove.related.FID_IID.txt \\
          /n/groups/price/UKBiobank/download_500K/w14048_20200204.FID_IID.txt \\
--exclude   /n/groups/price/UKBiobank/snpQC/autosome_maf_lt_0.001.txt \\
          /n/groups/price/UKBiobank/snpQC/autosome_missing_gt_0.1.txt \\
          /home/[username]/pleiotropy_o2/gwas/data/snps/height/coding.snps.for.height.tsv \\
--covar     /n/scratch3/users/n/[username]/gwas/data/snps/height/height.covariates.w.snps.chr$chr.tsv \\
--vif       1000000000 \\
--max-corr  1 \\
--pheno     /home/[username]/pleiotropy_o2/gwas/data/height.14048_ids.tab \\
--glm       hide-covar \\
--covar-variance-standardize \\
--threads   1 \\
--out       /n/scratch3/users/n/[username]/gwas/results/height/ukbb.height.assoc.chr${chr}

EOF
done

scp ~/tmp/run.height.gwas.chr*.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..1}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.height.gwas.chr${chr}.sh
done
'
