#!/bin/bash

### Lift the data

cat > ~/tmp/lift.zhang.pos.to.hg38.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=25G

mkdir -p /net/data/GTEx/topmed_europeans/zhang_reference

awk 'FNR > 1 {print "chr" $9 "\t" $10-1 "\t" $10 "\t" $11 "_" $12 "\tAF_" $13}' \
<(unzip -p /net/data/GTEx/gwas/zhang_2020_bc/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt.zip) \
> /net/data/GTEx/topmed_europeans/zhang_reference/zhang.hg19.pos.bed

/net/bin/x86_64/liftOver \
/net/data/GTEx/topmed_europeans/zhang_reference/zhang.hg19.pos.bed \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/topmed_europeans/zhang_reference/zhang.hg38.pos.bed \
/net/data/GTEx/topmed_europeans/zhang_reference/zhang.unlifted.hg19.to.38.bed

EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.zhang.pos.to.hg38.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/lift.zhang.pos.to.hg38.sh

ssh -J nconnally@ika.bwh.harvard.edu nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=20G /net/data/GTEx/topmed_europeans/tmp_src/lift.zhang.pos.to.hg38.sh"

### Perform the extraction

for chr in {1..22}; do
cat > ~/tmp/extract.zhang.snps.topmed.chr$chr.sh <<EOF
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0

mkdir -p /net/data/GTEx/topmed_europeans/zhang_hg38_plink

plink --vcf /net/data/GTEx/topmed_europeans/liu_hg38_bcf/topmed.eur.chr$chr.liu.vcf.gz \\
--allow-extra-chr \\
--chr $chr \\
--extract range <(awk '\$1 == "chr$chr" {print \$1 "\t" \$2+1 "\t" \$2+1 "\tinterval" NR}' /net/data/GTEx/topmed_europeans/zhang_reference/zhang.hg38.pos.bed) \\
--make-bed \\
--out /net/data/GTEx/topmed_europeans/zhang_hg38_plink/zhang.snps.topmed.chr$chr
EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/extract.zhang.snps.topmed.chr*.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/extract.zhang.snps.topmed.chr$chr.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'for chr in {1..22}; do
~/src/sgesh qsub /net/data/GTEx/topmed_europeans/tmp_src/extract.zhang.snps.topmed.chr$chr.sh
done'

### Redo the TOPMed names

cat > ~/tmp/redo.topmed.zhang.snp.names.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=10G
#$ -p 0

cd /net/data/GTEx/topmed_europeans/zhang_hg38_plink

for chr in {1..22}; do
plink --bfile zhang.snps.topmed.chr$chr \
--chr $chr \
--snps-only \
--maf 0.00001 \
--set-missing-var-ids @:#:\$1:\$2 \
--make-bed \
--out renamed.zhang.snps.chr$chr.topmed
done
EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/redo.topmed.zhang.snp.names.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/eo_files

rm ~/tmp/redo.topmed.zhang.snp.names.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'~/src/sgesh qsub /net/data/GTEx/topmed_europeans/eo_files/redo.topmed.zhang.snp.names.sh'

### Liftover coding SNPs

cat > ~/tmp/lift.zhang.coding.snps.to.hg38.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=25G

/net/bin/x86_64/liftOver \
<(awk -F ':' '{print "chr" $1 "\t" $2-1 "\t" $2 "\t" $3 ":" $4}' /net/data/GTEx/gwas/zhang_2020_bc/cojo/coding_snps/coding.snps.for.bc.tsv) \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/coding_snps/coding.snps.for.bc.hg38.bed \
/net/data/GTEx/topmed_europeans/zhang_reference/coding.snps.for.bc.unlifted.hg19.to.38.bed

awk '{gsub(/chr/, ""); print $1 ":" $2+1 ":" $4}' \
/net/data/GTEx/gwas/zhang_2020_bc/cojo/coding_snps/coding.snps.for.bc.hg38.bed \
> /net/data/GTEx/gwas/zhang_2020_bc/cojo/coding_snps/coding.snps.for.bc.hg38.tsv

EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.zhang.coding.snps.to.hg38.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/lift.zhang.coding.snps.to.hg38.sh

ssh -J nconnally@ika.bwh.harvard.edu nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=20G /net/data/GTEx/topmed_europeans/tmp_src/lift.zhang.coding.snps.to.hg38.sh"
