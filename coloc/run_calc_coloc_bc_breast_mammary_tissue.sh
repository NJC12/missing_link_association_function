#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/bc/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.bc.breast_mammary_tissue.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.R \\
gwas/zhang_2020_bc/cojo/results_hg19/regions/zhang.bc.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/breast_mammary_tissue/regions/breast_mammary_tissue.$chr.$start.$end.txt \\
coloc/bc/regions/coloc.bc.breast_mammary_tissue.$chr.$start.$end.txt \\
breast_mammary_tissue \\
247173

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.bc.breast_mammary_tissue.sh

done < gwas/zhang_2020_bc/cojo/results_hg19/main_pheno_peaks/zhang_bc_hg19.indexSNP.tsv
