#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/hdl/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.hdl.liver.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.ukbb.R \\
gwas/ukbb_hdl_conditioned/regions/ukbb.hdl.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/liver/regions/liver.$chr.$start.$end.txt \\
coloc/hdl/regions/coloc.hdl.liver.$chr.$start.$end.txt \\
liver \\
337489

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.hdl.liver.sh

done < gwas/ukbb_hdl_conditioned/main_pheno_peaks_hg19/ukbb_hdl_conditioned_hg19.indexSNP.tsv
