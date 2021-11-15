#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/cd/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.cd.colon_transverse.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.R \\
gwas/liu_2015_and_goyette_2015_ibd/cojo/results_hg19/regions/liu.cd.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/colon_transverse/regions/colon_transverse.$chr.$start.$end.txt \\
coloc/cd/regions/coloc.cd.colon_transverse.$chr.$start.$end.txt \\
colon_transverse \\
86640  

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.cd.colon_transverse.sh

done < gwas/liu_2015_and_goyette_2015_ibd/cojo/results_hg19/main_pheno_peaks/liu_cd_hg19.indexSNP.tsv
