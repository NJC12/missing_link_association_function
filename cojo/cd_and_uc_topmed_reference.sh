#!/bin/bash

### Extract the SNPs used in these GWAS from the 5,000 TOPMed Europeans we are using

for chr in {1..22}; do
cat > ~/tmp/tmp.extract.liu.snps.topmed.chr$chr.sh <<EOF
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=20G,h_vmem=20G
#$ -p 0



mkdir -p /net/data/GTEx/topmed_europeans/liu_hg38_plink

plink --vcf /net/data/GTEx/topmed_europeans/liu_hg38_bcf/topmed.eur.chr$chr.liu.vcf.gz \\
--allow-extra-chr \\
--chr $chr \\
--extract range <(awk '\$1 == "chr$chr" {print \$1 "\t" \$2+1 "\t" \$2+1 "\tinterval" NR}' /net/data/GTEx/topmed_europeans/liu_reference/liu.hg38.pos.bed) \\
--make-bed \\
--out /net/data/GTEx/topmed_europeans/liu_hg38_plink/liu.snps.topmed.chr$chr
EOF

done

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/tmp.extract.liu.snps.topmed.chr*.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/tmp.extract.liu.snps.topmed.chr$chr.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'for chr in {1..22}; do
~/src/sgesh qsub /net/data/GTEx/topmed_europeans/tmp_src/tmp.extract.liu.snps.topmed.chr$chr.sh
done'

### Change the TOPMed SNP names to match


cat > ~/tmp/redo.topmed.snp.names.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/gwas/eo_files
#$ -o /net/data/GTEx/gwas/eo_files
#$ -l mf=10G
#$ -p 0

cd /net/data/GTEx/topmed_europeans/liu_hg38_plink

for chr in {1..22}; do
plink --bfile liu.snps.topmed.chr$chr \
--chr $chr \
--snps-only \
--maf 0.00001 \
--set-missing-var-ids @:#:\$1:\$2 \
--make-bed \
--out renamed.liu.snps.chr$chr.topmed
done
EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/redo.topmed.snp.names.sh \
nconnally@tako.partners.org:/net/data/GTEx/gwas/tmp_src

rm ~/tmp/redo.topmed.snp.names.sh

ssh -J nconnally@ika.bwh.harvard.edu \
nconnally@tako.partners.org \
'~/src/sgesh qsub /net/data/GTEx/gwas/tmp_src/redo.topmed.snp.names.sh'

### Liftover coding SNPs

cat > ~/tmp/lift.liu.coding.snps.to.hg38.sh <<'EOF'
#!/bin/bash
#$ -e /net/data/GTEx/topmed_europeans/eo_files
#$ -o /net/data/GTEx/topmed_europeans/eo_files
#$ -l mf=25G



/net/bin/x86_64/liftOver \
<(awk -F ':' '{print "chr" $1 "\t" $2-1 "\t" $2 "\t" $3 ":" $4}' /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/coding_snps/coding.snps.for.ibd.tsv) \
/net/home/nconnally/data/liftover_chains/hg19ToHg38.over.chain.gz \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/coding_snps/coding.snps.for.ibd.hg38.bed \
/net/data/GTEx/topmed_europeans/liu_reference/coding.snps.for.ibd.unlifted.hg19.to.38.bed

awk '{gsub(/chr/, ""); print $1 ":" $2+1 ":" $4}' \
/net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/coding_snps/coding.snps.for.ibd.hg38.bed \
> /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/coding_snps/coding.snps.for.ibd.hg38.tsv

EOF

scp -o ProxyJump=nconnally@ika.bwh.harvard.edu \
~/tmp/lift.liu.coding.snps.to.hg38.sh \
nconnally@tako.partners.org:/net/data/GTEx/topmed_europeans/tmp_src

rm ~/tmp/lift.liu.coding.snps.to.hg38.sh

ssh -J nconnally@ika.bwh.harvard.edu nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=20G /net/data/GTEx/topmed_europeans/tmp_src/lift.liu.coding.snps.to.hg38.sh"
