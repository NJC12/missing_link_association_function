#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/height/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.height.muscle_skeletal.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.ukbb.R \\
gwas/ukbb_height_conditioned/regions/ukbb.height.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/muscle_skeletal/regions/muscle_skeletal.$chr.$start.$end.txt \\
coloc/height/regions/coloc.height.muscle_skeletal.$chr.$start.$end.txt \\
muscle_skeletal \\
337489

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.height.muscle_skeletal.sh

done < gwas/ukbb_height_conditioned/main_pheno_peaks_hg19/ukbb_height_conditioned_hg19.indexSNP.tsv
