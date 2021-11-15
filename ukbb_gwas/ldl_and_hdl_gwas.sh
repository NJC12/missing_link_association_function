#!/bin/bash

### Getting the data in the correct formats

cat > ~/tmp/create.cholesterol.pheno.file.sh <<'EOF'
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/create.cholesterol.pheno.file.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/create.cholesterol.pheno.file.sh.o%j
#SBATCH --mem=15G
#SBATCH -p short
#SBATCH -t 5



grep -f <(grep "1140861922" /n/groups/price/UKBiobank/app10438assoc/ukb4986.tab | cut -f 1) \
/n/groups/price/UKBiobank/app10438assoc/ukb4989.IDmap.tab | \
cut -f 2 \
> /home/[username]/pleiotropy_o2/gwas/data/subjects.on.lipid.lowering.medication.txt

awk '
NR == FNR && NR > 1 {
key[$1]=$2
}
NR != FNR && FNR == 1 {
for (i=1; i<=NF; i++) {
  if ($i == "f.30760.0.0") {c1=i};
  if ($i == "f.30780.0.0") {c2=i};
}
print "FID\tIID\t" $c1 "\t" $c2
}
NR != FNR && FNR > 1 {
print key[$1] "\t" key[$1] "\t" $c1 "\t" $c2
}' \
/n/groups/price/UKBiobank/app10438assoc/ukb4989.IDmap.tab \
/n/groups/price/UKBiobank/app10438assoc/ukb27621.tab | \
grep -vf /home/[username]/pleiotropy_o2/gwas/data/subjects.on.lipid.lowering.medication.txt \
> /home/[username]/pleiotropy_o2/gwas/data/cholesterol.14048_ids.tab

# awk '
# FNR == NR {meds[$1]=1}
# FNR != NR {print $1 "\t" $2 "\t" $5 "\t" $20 "\t" $36 "\t" $73 "\t" $74 "\t" $75 "\t" $76 "\t" $77 "\t" $78 "\t" $79 "\t" $80 "\t" $81 "\t" $82 "\t" 0+meds[$1]}
# ' \
# /home/[username]/pleiotropy_o2/gwas/data/subjects.on.lipid.lowering.medication.txt \
# <(zcat /n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz) \
# > /home/[username]/pleiotropy_o2/gwas/data/cholesterol.covariates.14048_ids.tab

awk '
{print $1 "\t" $2 "\t" $5 "\t" $20 "\t" $21 "\t" $36 "\t" $73 "\t" $74 "\t" $75 "\t" $76 "\t" $77 "\t" $78 "\t" $79 "\t" $80 "\t" $81 "\t" $82}
' \
<(zcat /n/groups/price/UKBiobank/app10438assoc/ukb4777.processed_and_post.plinkPCs.tab.gz) \
> /home/[username]/pleiotropy_o2/gwas/data/cholesterol.covariates.14048_ids.tab
EOF

scp ~/tmp/create.cholesterol.pheno.file.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/create.cholesterol.pheno.file.sh
'

### Extract cholesterol coding SNPs

mem_array=( 1000 256 256 256 256 256 256 256 256 150 150 150 150 150 150 150 150 100 150 100 100 100 100 )
time_array=( 1000 6 6 6 6 6 6 6 6 4 4 4 4 4 4 4 4 2 6 2 2 2 2 )

for chr in {1..22}; do

cat > ~/tmp/create.cholesterol.pheno.file.chr$chr.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/create.cholesterol.pheno.file.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/create.cholesterol.pheno.file.sh.o%j
#SBATCH --mem=${mem_array[$chr]}G
#SBATCH -p short
#SBATCH -t ${time_array[$chr]}:00:00



mkdir -p /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol

module load plink2/2.0

plink2 \\
--bgen       /n/groups/price/UKBiobank/bgen_MAF001_500K_v3/UKB_MAF0.001_v3.$chr.bgen \\
--sample     /n/groups/price/UKBiobank/download_500K/ukb14048_imp_chr1_v3_s487395.sample \\
--remove     /n/groups/price/UKBiobank/sampleQC/remove.nonStringentBritish.FID_IID.txt \\
             /n/groups/price/UKBiobank/sampleQC/remove.related.FID_IID.txt \\
             /n/groups/price/UKBiobank/download_500K/w14048_20200204.FID_IID.txt \\
--extract    /home/[username]/pleiotropy_o2/gwas/data/snps/cholesterol/coding.snps.for.cholesterol.tsv \\
--export     vcf \\
--out        /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/ukbb.cholesterol.codings.snps.chr$chr

# mv /n/scratch3/users/n/[username]/data/snps/cholesterol/ukbb.cholesterol.codings.snps.chr$chr.vcf \
# /home/[username]/pleiotropy_o2/gwas/data/snps/cholesterol/

EOF

done

scp ~/tmp/create.cholesterol.pheno.file.chr*.sh [username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..22}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/create.cholesterol.pheno.file.chr$chr.sh
done
'

for chr in {1..22}; do
cat > ~/tmp/run.transpose.vcfs.chr$chr.py <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.transpose.cholesterol.vcfs.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.transpose.cholesterol.vcfs.sh.o%j
#SBATCH --mem=100G
#SBATCH -p short
#SBATCH -t 01:00:00

mkdir -p /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol

module load gcc/6.2.0
module load python/3.7.4

python3.7 /home/[username]/pleiotropy_o2/gwas/tmp_src/transpose.vcfs.py \
/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/ukbb.cholesterol.codings.snps.chr$chr.vcf \
/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf

EOF
done

scp \
~/tmp/transpose.vcfs.py \
~/tmp/run.transpose.vcfs.chr*.py \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..22}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.transpose.vcfs.chr$chr.py
done
'

### Remove correlated covariates

cat > ~/tmp/remove.correlated.covariates.R <<'EOF'


args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]

chr <- read.table(infile, header=F, as.is=T, sep='\t')
names <- chr[, 1:2]
snps <- chr[, 3:ncol(chr)]

cors <- cor(snps)
diag(cors) <- 0
round <- 0
while(any(cors > 0.9 | cors < -0.9, na.rm=T)) {
  round <- round + 1
  print(round)
  high.cors <- apply(cors, 2, function(x) sum(x > 0.9 | x < -0.9, na.rm=T))
  snps <- snps[, -which.max(high.cors)]
  print(paste('Dim:', dim(snps)))
  cors <- cor(snps)
  diag(cors) <- 0
}

write.table(cbind(names, snps), outfile, row.names=F, col.names=F, quote=F, sep='\t')

EOF

chr=6
cat > ~/tmp/run.remove.correlated.covariates.chr$chr.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.remove.correlated.covariates.chr$chr.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.remove.correlated.covariates.chr$chr.sh.o%j
#SBATCH --mem=20G
#SBATCH -p short
#SBATCH -t 12:00:00



module load gcc/6.2.0
module load R/3.6.1

# file=/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/BEFORE.COR.REMOVED.transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf
# [ -f "$file" ] && echo "File without correlation removed already exists" || \\
# sed -r 's/_/\t/g' /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf \\
# > /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/BEFORE.COR.REMOVED.transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf

Rscript /home/[username]/pleiotropy_o2/gwas/tmp_src/remove.correlated.covariates.R \\
/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/BEFORE.COR.REMOVED.transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf \\
/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf
EOF


scp \
~/tmp/remove.correlated.covariates.R \
~/tmp/run.remove.correlated.covariates.chr$chr.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
chr=6
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.remove.correlated.covariates.chr$chr.sh
'

### Combine the regular covariates with the coding SNPs

for chr in {1..22}; do
cat > ~/tmp/run.combine.cholesterol.covariates.chr$chr.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.combine.cholesterol.covariates.chr$chr.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.combine.cholesterol.covariates.chr$chr.sh.o%j
#SBATCH --mem=5G
#SBATCH -p short
#SBATCH -t 00:10:00

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
/home/[username]/pleiotropy_o2/gwas/data/cholesterol.covariates.14048_ids.tab \\
/n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/transposed.ukbb.cholesterol.codings.snps.chr$chr.vcf | \\
sort | uniq \\
> /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/cholesterol.covariates.w.snps.chr$chr.tsv

EOF
done

scp \
~/tmp/run.combine.cholesterol.covariates.chr*.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..22}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.combine.cholesterol.covariates.chr$chr.sh
done
'

### Run the GWAS

for chr in {1..22}; do
cat > ~/tmp/run.cholesterol.gwas.chr${chr}.sh <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -e /home/[username]/pleiotropy_o2/gwas/eo_files/run.cholesterol.gwas.chr${chr}.sh.e%j
#SBATCH -o /home/[username]/pleiotropy_o2/gwas/eo_files/run.cholesterol.gwas.chr${chr}.sh.o%j
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=100G
#SBATCH -p medium
#SBATCH -t 5-00:00:00



mkdir -p /n/scratch3/users/n/[username]/gwas/results/cholesterol

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
--covar     /n/scratch3/users/n/[username]/gwas/data/snps/cholesterol/cholesterol.covariates.w.snps.chr$chr.tsv \\
--vif       1000000000 \\
--max-corr  1 \\
--pheno     /home/[username]/pleiotropy_o2/gwas/data/cholesterol.14048_ids.tab \\
--glm       hide-covar \\
--covar-variance-standardize \\
--threads   1 \\
--out       /n/scratch3/users/n/[username]/gwas/results/cholesterol/ukbb.cholesterol.assoc.chr${chr}

EOF
done

scp ~/tmp/run.cholesterol.gwas.chr*.sh \
[username]@[cluster]:/home/[username]/pleiotropy_o2/gwas/tmp_src

ssh [username]@[cluster] '
for chr in {1..22}; do
sbatch /home/[username]/pleiotropy_o2/gwas/tmp_src/run.cholesterol.gwas.chr${chr}.sh
done
'
