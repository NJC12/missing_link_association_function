#$ -S /bin/bash
#$ -e /net/data/GTEx/eo_files
#$ -o /net/data/GTEx/eo_files
#$ -l mf=15G
#$ -V

cd /net/data/GTEx/gwas/ukbb_ldl_conditioned

mkdir -p regions

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

awk -v chr=$chr -v start=$start -v end=$end \
'NR == 1 || ($1 == chr && $2 >= start && $2 <= end)' \
ukbb.ldl.conditioned.hg19.renamed.tsv \
> regions/ukbb.ldl.$chr.$start.$end.tsv

done < main_pheno_peaks_hg19/ukbb_ldl_conditioned_hg19.indexSNP.tsv
