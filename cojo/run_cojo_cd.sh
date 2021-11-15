#!/bin/bash

### Creating .ma reference files for COJO.

cat > ~/tmp/prepare.liu.cd.for.cojo.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/eo_files
#$ -o /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/eo_files
#$ -l mf=5G,h_vmem=5G
#$ -p 0

for chr in {1..22}; do

awk -v chr=$chr 'NR > 1 && $1 == chr {print "chr" $1 "\t" $3-1 "\t" $3 "\t" $4 "_" $5 "_" $7 "_" log($9) "_" $10 "_" $11 "_20883";
print "chr" $1 "\t" $3-1 "\t" $3 "\t" $5 "_" $4 "_" 1-$7 "_" log($9) "_" $10 "_" $11 "_20883"}' \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/EUR.CD.gwas_info03_filtered.assoc \
> /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.hg19.chr$chr.bed

/net/bin/x86_64/liftOver \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.hg19.chr$chr.bed \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.hg38.chr$chr.bed \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.unlifted.hg19.to.38.chr$chr.bed

awk -v chr=$chr '$1 == "chr"chr' \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.hg38.chr$chr.bed | \
awk -F '[\t_]' 'BEGIN {print "SNP\tA1\tA2\tfreq\tb\tse\tp\tN"} {gsub(/chr/, ""); print $1 ":" $2+1 ":" $4 ":" $5 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' \
> /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/ma_files/liu.cd.hg38.chr$chr.ma

done
EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/prepare.liu.cd.for.cojo.sh \
nconnally@tako.partners.org:/net/data/GTEx/gwas/tmp_src

rm ~/tmp/prepare.liu.cd.for.cojo.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'~/src/sgesh qsub /net/data/GTEx/gwas/tmp_src/prepare.liu.cd.for.cojo.sh'

### Running COJO

dir="/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo"
study="liu"
study2=".cd"
pheno="ibd"

for chr in {1..22}; do

cat > ~/tmp/run.cojo.liu.cd.chr$chr.sh <<EOF
#!/bin/bash
#$ -e $dir/eo_files
#$ -o $dir/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0

gcta64 --bfile /net/data/GTEx/topmed_europeans/${study}_hg38_plink/renamed.$study.snps.chr$chr.topmed \\
--cojo-file $dir/ma_files/$study$study2.hg38.chr$chr.ma \\
--extract <(grep "^${chr}:" $dir/coding_snps/coding.snps.for.$pheno.hg38.tsv) \\
--cojo-slct \\
--cojo-p 5e-2 \\
--out $dir/sig_coding_snps/$study$study2.chr$chr.cond.sig.coding.snps

gcta64 --bfile /net/data/GTEx/topmed_europeans/${study}_hg38_plink/renamed.$study.snps.chr$chr.topmed \\
--cojo-file $dir/ma_files/$study$study2.hg38.chr$chr.ma \\
--cojo-cond <(awk '\$8 <= 0.05 {print \$2}' $dir/sig_coding_snps/$study$study2.chr$chr.cond.sig.coding.snps.jma.cojo) \\
--out $dir/results/$study$study2.chr$chr.output.conditioned

EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/run.cojo.liu.cd.chr*.sh \
nconnally@tako.partners.org:$dir/tmp_src

rm ~/tmp/run.cojo.liu.cd.chr*.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
"for chr in {1..22}; do
~/src/sgesh qsub $dir/tmp_src/run.cojo.liu.cd.chr\$chr.sh
done"

### Lift results back to hg19 for our study

dir="/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo"
study="liu"
pheno="cd"

for chr in {1..22}; do

cat > ~/tmp/lift.liu.to.hg19.sh <<EOF
#!/bin/bash
#$ -e $dir/eo_files
#$ -o $dir/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0

mkdir -p $dir/results_hg19

printf "" > $dir/results/$study.$pheno.conditioned.tsv
for chr in {1..22}; do awk 'NR > 1' $dir/results/$study.$pheno.chr\$chr.output.conditioned.cma.cojo >> $dir/results/$study.$pheno.conditioned.tsv; done

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
~/tmp/lift.liu.to.hg19.sh \
nconnally@tako.partners.org:$dir/tmp_src

rm ~/tmp/lift.liu.to.hg19.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
"~/src/sgesh qsub $dir/tmp_src/lift.liu.to.hg19.sh"
