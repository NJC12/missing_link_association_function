#!/bin/bash

### Extract the correct SNPs from TOPMed

cat > ~/tmp/lift.mahajan.pos.to.hg38.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=25G



mkdir -p /net/data/GTEx/topmed_europeans/mahajan_reference

awk 'FNR > 1 {print "chr" $2 "\t" $3-1 "\t" $3 "\t" $4 "_" $5 "\t" $6}' \
/net/data/GTEx/gwas/mahajan_2018_t2d/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt \
> /net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.hg19.pos.bed

/net/bin/x86_64/liftOver \
/net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.hg19.pos.bed \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.hg38.pos.bed \
/net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.unlifted.hg19.to.38.bed

EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.mahajan.pos.to.hg38.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/lift.mahajan.pos.to.hg38.sh

ssh -J nconnally@ika.bwh.harvard.edu nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=20G /net/data/GTEx/topmed_europeans/tmp_src/lift.mahajan.pos.to.hg38.sh"

for chr in {1..22}; do
cat > ~/tmp/extract.mahajan.snps.topmed.chr$chr.sh <<EOF
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0



mkdir -p /net/data/GTEx/topmed_europeans/mahajan_hg38_plink

plink --vcf /net/data/GTEx/topmed_europeans/liu_hg38_bcf/topmed.eur.chr$chr.liu.vcf.gz \\
--allow-extra-chr \\
--chr $chr \\
--extract range <(awk '\$1 == "chr$chr" {print \$1 "\t" \$2+1 "\t" \$2+1 "\tinterval" NR}' /net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.hg38.pos.bed) \\
--make-bed \\
--out /net/data/GTEx/topmed_europeans/mahajan_hg38_plink/mahajan.snps.topmed.chr$chr
EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/extract.mahajan.snps.topmed.chr*.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/extract.mahajan.snps.topmed.chr$chr.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'for chr in {1..22}; do
~/src/sgesh qsub /net/data/GTEx/topmed_europeans/tmp_src/extract.mahajan.snps.topmed.chr$chr.sh
done'

### Liftover coding SNPs

cat > ~/tmp/lift.mahajan.coding.snps.to.hg38.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=25G



mkdir -p /net/data/GTEx/topmed_europeans/mahajan_reference

awk 'FNR > 1 {print "chr" $2 "\t" $3-1 "\t" $3 "\t" $4 "_" $5 "\t" $6}' \
/net/data/GTEx/gwas/mahajan_2018_t2d/Mahajan.NatGenet2018b.T2Dbmiadj.European.txt \
> /net/data/GTEx/topmed_europeans/mahajan_reference/mahajan.hg19.pos.bed

/net/bin/x86_64/liftOver \
<(awk -F ':' '{print "chr" $1 "\t" $2-1 "\t" $2 "\t" $3 ":" $4}' /net/data/GTEx/gwas/mahajan_2018_t2d/cojo/coding_snps/coding.snps.for.t2d.tsv) \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/gwas/mahajan_2018_t2d/cojo/coding_snps/coding.snps.for.t2d.hg38.bed \
/net/data/GTEx/topmed_europeans/mahajan_reference/coding.snps.for.t2d.unlifted.hg19.to.38.bed

awk '{gsub(/chr/, ""); print $1 ":" $2+1 ":" $4}' \
/net/data/GTEx/gwas/mahajan_2018_t2d/cojo/coding_snps/coding.snps.for.t2d.hg38.bed \
> /net/data/GTEx/gwas/mahajan_2018_t2d/cojo/coding_snps/coding.snps.for.t2d.hg38.tsv

EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.mahajan.coding.snps.to.hg38.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/lift.mahajan.coding.snps.to.hg38.sh

ssh -J nconnally@ika.bwh.harvard.edu nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=20G /net/data/GTEx/topmed_europeans/tmp_src/lift.mahajan.coding.snps.to.hg38.sh"
