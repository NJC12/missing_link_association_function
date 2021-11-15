######################################
### Extract 5,000 TOPMed Europeans ###
######################################

# Extracts all SNPs, because each GWAS has a slightly different SNP set.

for chr in {1..22}; do
cat > ~/tmp/extract.chr$chr.eurs.all.snps.vcf.sh <<EOF
#!/bin/bash
#$ -e /net/home/nconnally/data/topmed_europeans/eo_files
#$ -o /net/home/nconnally/data/topmed_europeans/eo_files



PATH=\$PATH:/net/bin/x86_64/

mkdir -p /net/home/nconnally/data/topmed_europeans/hg38_all_snps
cd /net/home/nconnally/data/topmed_europeans/hg38_all_snps

bcftools view /net/data/pub/topmed/TOPMed_Freeze_6_Phased/freeze.6a.chr$chr.pass_only.phased.bcf \\
-S /net/home/nconnally/data/topmed_europeans/reference/topmed.4700.europeans.list.txt \\
-O v \\
> topmed.hg38.eurs.all.snps.chr$chr.vcf
EOF

scp ~/tmp/extract.chr$chr.eurs.all.snps.vcf.sh \
nconnally@tako.partners.org:/net/home/nconnally/data/topmed_europeans/tmp_src
rm ~/tmp/extract.chr$chr.eurs.all.snps.vcf.sh

done

for chr in {1..22}; do
ssh nconnally@tako.partners.org \
"~/src/sgesh qsub -p 0 -l mf=10G /net/home/nconnally/data/topmed_europeans/tmp_src/extract.chr$chr.eurs.all.snps.vcf.sh"

done
