#$ -S /bin/bash
#$ -e /net/data/GTEx/eo_files
#$ -o /net/data/GTEx/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx/gwas/zhang_2020_bc/cojo/results_hg19

mkdir -p regions

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

awk -v chr=$chr -v start=$start -v end=$end \
'NR == 1 || ($1 == chr && $3 >= start && $3 <= end)' \
zhang.bc.hg19.conditioned.cojo \
> regions/zhang.bc.$chr.$start.$end.tsv

done < main_pheno_peaks/zhang_bc_hg19.indexSNP.tsv
