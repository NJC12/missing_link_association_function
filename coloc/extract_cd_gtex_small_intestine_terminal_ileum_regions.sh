#$ -S /bin/bash
#$ -e /net/data/GTEx/eo_files
#$ -o /net/data/GTEx/eo_files
#$ -l mf=15G
#$ -V

mkdir -p /net/data/GTEx/GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/small_intestine_terminal_ileum/regions
cd /net/data/GTEx/GTEx_Analysis_v7_QTLs/GTEx_Analysis_v7_eQTL_all_associations/small_intestine_terminal_ileum

while read -r line; do

chr=$(echo $line | awk '{print $1}')
if [ $chr = "CHR" ]; then
continue
fi
start=$(echo $line | awk '{print $4}')
end=$(echo $line | awk '{print $5}')

awk -v chr=$chr -v start=$start -v end=$end '\
NR == 1 {print $0}
NR > 1 {split($2,pos,"_");
if (pos[1] == chr && pos[2] >= start && pos[2] <= end) {print $0}}' \
Small_Intestine_Terminal_Ileum.chr$chr.txt \
> regions/small_intestine_terminal_ileum.$chr.$start.$end.txt

done < /net/data/GTEx/gwas/liu_2015_and_goyette_2015_ibd/cojo/results_hg19/main_pheno_peaks/liu_cd_hg19.indexSNP.tsv
