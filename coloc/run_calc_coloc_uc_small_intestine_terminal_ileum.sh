#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx

mkdir -p coloc/uc/regions/

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

cat > /net/data/GTEx/src/calc.coloc.$chr.$start.$end.uc.small_intestine_terminal_ileum.sh <<EOF2
#$ -S /bin/bash
#$ -e /net/data/GTEx/coloc/eo_files
#$ -o /net/data/GTEx/coloc/eo_files
#$ -l mf=7G
#$ -V

cd /net/data/GTEx

Rscript /net/data/GTEx/src/calc.coloc.piecemeal.R \\
gwas/liu_2015_and_goyette_2015_ibd/cojo/results_hg19/regions/liu.uc.$chr.$start.$end.tsv \\
GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/small_intestine_terminal_ileum/regions/small_intestine_terminal_ileum.$chr.$start.$end.txt \\
coloc/uc/regions/coloc.uc.small_intestine_terminal_ileum.$chr.$start.$end.txt \\
small_intestine_terminal_ileum \\
86640  

EOF2

qsub -p 0 -l mf=7G /net/data/GTEx/src/calc.coloc.$chr.$start.$end.uc.small_intestine_terminal_ileum.sh

done < gwas/liu_2015_and_goyette_2015_ibd/cojo/results_hg19/main_pheno_peaks/liu_uc_hg19.indexSNP.tsv
