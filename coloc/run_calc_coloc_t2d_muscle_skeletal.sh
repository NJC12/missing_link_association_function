#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/t2d/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.t2d.muscle_skeletal.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.R \\
gwas/mahajan_2018_t2d/cojo/results_hg19/regions/mahajan.t2d.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/muscle_skeletal/regions/muscle_skeletal.$chr.$start.$end.txt \\
coloc/t2d/regions/coloc.t2d.muscle_skeletal.$chr.$start.$end.txt \\
muscle_skeletal \\
132767

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.t2d.muscle_skeletal.sh

done < gwas/mahajan_2018_t2d/cojo/results_hg19/main_pheno_peaks/mahajan_t2d_hg19.indexSNP.tsv
