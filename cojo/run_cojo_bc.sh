#!/bin/bash

### Create the .ma reference needed for COJO

cat > ~/tmp/prepare.zhang.for.cojo.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/gwas/zhang_2020_bc/cojo/eo_files
#$ -o /net/data/GTEx/gwas/zhang_2020_bc/cojo/eo_files
#$ -l mf=5G,h_vmem=5G
#$ -p 0

for chr in {1..22}; do

awk -v chr=$chr 'NR > 1 && $9 == chr {
print "chr" $9 "\t" $10-1 "\t" $10 "\t" $40 "_" $41 "_"   $29 "_" $42 "_" $44 "_" $46 "_247173";
print "chr" $9 "\t" $10-1 "\t" $10 "\t" $41 "_" $40 "_" 1-$29 "_" $42 "_" $44 "_" $46 "_247173"}' \
<(unzip -p /net/data/GTEx/gwas/zhang_2020_bc/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt.zip) \
> /net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.hg19.chr$chr.bed

/net/bin/x86_64/liftOver \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.hg19.chr$chr.bed \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.hg38.chr$chr.bed \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.unlifted.hg19.to.38.chr$chr.bed

awk -v chr=$chr '$1 == "chr"chr' \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.hg38.chr$chr.bed | \
awk -F '[\t_]' 'BEGIN {print "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"} {gsub(/chr/, ""); print $1 ":" $2+1 ":" $4 ":" $5 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
> /net/data/GTEx/gwas/zhang_2020_bc/cojo/ma_files/zhang.hg38.chr$chr.ma

done
EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/prepare.zhang.for.cojo.sh \
nconnally@tako.partners.org:/net/data/GTEx/gwas/tmp_src

rm ~/tmp/prepare.zhang.for.cojo.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'~/src/sgesh qsub /net/data/GTEx/gwas/tmp_src/prepare.zhang.for.cojo.sh'

### Run COJO

dir="/net/data/GTEx/gwas/zhang_2020_bc/cojo"
study="zhang"
pheno="bc"

for chr in {1..22}; do

cat > ~/tmp/run.cojo.zhang.chr$chr.sh <<EOF
#!/bin/bash
#$ -e $dir/eo_files
#$ -o $dir/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0

gcta64 --bfile /net/data/GTEx/topmed_europeans/${study}_hg38_plink/renamed.$study.snps.chr$chr.topmed \\
--cojo-file $dir/ma_files/$study.hg38.chr$chr.ma \\
--extract <(grep "^${chr}:" $dir/coding_snps/coding.snps.for.$pheno.hg38.tsv) \\
--cojo-slct \\
--cojo-p 5e-2 \\
--out $dir/sig_coding_snps/$study.chr$chr.cond.sig.coding.snps

gcta64 --bfile /net/data/GTEx/topmed_europeans/${study}_hg38_plink/renamed.$study.snps.chr$chr.topmed \\
--cojo-file $dir/ma_files/$study.hg38.chr$chr.ma \\
--cojo-cond <(awk '\$8 <= 0.05 {print \$2}' $dir/sig_coding_snps/$study.chr$chr.cond.sig.coding.snps.jma.cojo) \\
--out $dir/results/$study.chr$chr.output.conditioned

EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/run.cojo.zhang.chr*.sh \
nconnally@tako.partners.org:$dir/tmp_src

rm ~/tmp/run.cojo.zhang.chr*.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
"for chr in {1..22}; do
~/src/sgesh qsub $dir/tmp_src/run.cojo.zhang.chr\$chr.sh
done"

### Lift back to hg19 for our study

dir="/net/data/GTEx/gwas/zhang_2020_bc/cojo"
study="zhang"
pheno="bc"

for chr in {1..22}; do

cat > ~/tmp/lift.zhang.to.hg19.sh <<EOF
#!/bin/bash
#$ -e $dir/eo_files
#$ -o $dir/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0

mkdir -p $dir/results_hg19

printf "" > $dir/results/$study.$pheno.conditioned.tsv
for chr in {1..22}; do awk 'NR > 1' $dir/results/$study.chr\$chr.output.conditioned.cma.cojo >> $dir/results/$study.$pheno.conditioned.tsv; done

awk -F '[\t:]' 'NR > 1 {print "chr"\$1 "\t" \$6-1 "\t" \$6 "\t" \$4 "_" \$5 "_" \$7 "_" \$8 "_" \$9 "_" \$10 "_" \$11 "_" \$12 "_" \$13 "_" \$14 "_" \$15 "_" \$16}' \
$dir/results/$study.$pheno.conditioned.tsv \
> $dir/results_hg19/$study.$pheno.hg38.conditioned.bed

/net/bin/x86_64/liftOver \
$dir/results_hg19/$study.$pheno.hg38.conditioned.bed \
/net/home/nconnally/data/liftover_chains/hg38ToHg19.over.chain.gz \
$dir/results_hg19/$study.$pheno.hg19.conditioned.bed \
$dir/results_hg19/$study.unlifted.hg38.to.19.bed

awk -F '[\t_]' 'BEGIN {print "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp\tn\tfreq_geno\tbC\tbC_se\tpC"}
{gsub(/chr/, ""); print \$1 "\t" \$1 ":" \$2+1 ":" \$4 ":" \$5 "\t" \$2+1 "\t" \$6 "\t" \$7 "\t" \$8 "\t" \$9 "\t" \$10 "\t" \$11 "\t" \$12 "\t" \$13 "\t" \$14 "\t" \$15}' \
$dir/results_hg19/$study.$pheno.hg19.conditioned.bed \
> $dir/results_hg19/$study.$pheno.hg19.conditioned.cojo

EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.zhang.to.hg19.sh \
nconnally@tako.partners.org:$dir/tmp_src

rm ~/tmp/lift.zhang.to.hg19.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
"~/src/sgesh qsub $dir/tmp_src/lift.zhang.to.hg19.sh"
