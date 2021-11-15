#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/ldl/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.ldl.liver.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.ukbb.R \\
gwas/ukbb_ldl_conditioned/regions/ukbb.ldl.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/liver/regions/liver.$chr.$start.$end.txt \\
coloc/ldl/regions/coloc.ldl.liver.$chr.$start.$end.txt \\
liver \\
337489

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.ldl.liver.sh

done < gwas/ukbb_ldl_conditioned/main_pheno_peaks_hg19/ukbb_ldl_conditioned_hg19.indexSNP.tsv
