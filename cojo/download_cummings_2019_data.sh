#!/bin/bash

# Download the data from Cummings et al. 2019 from Google Cloud
gsutil -m cp gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.GTEx.v7.021520.tsv.bgz .
cp all.possible.snvs.tx_annotated.GTEx.v7.021520.tsv.bgz all.possible.snvs.tx_annotated.GTEx.v7.021520.tsv.gz
gunzip all.possible.snvs.tx_annotated.GTEx.v7.021520.tsv.gz

mkdir -p tmp_src eo_files

cat > ~/tmp/cummings.gtex.tissue.key.sh <<'EOF'
s/\{//g
s/\[//g
s/\}//g
s/\]//g
s/AdrenalGland/Adrenal_Gland/g
s/Brain_Anteriorcingulatecortex_BA24_/Brain_Anterior_cingulate_cortex_BA24_/g
s/Brain_Caudate_basalganglia_/Brain_Caudate_basal_ganglia_/g
s/Brain_CerebellarHemisphere/Brain_Cerebellar_Hemisphere/g
s/Brain_FrontalCortex_BA9_/Brain_Frontal_Cortex_BA9_/g
s/Brain_Nucleusaccumbens_basalganglia_/Brain_Nucleus_accumbens_basal_ganglia_/g
s/Brain_Putamen_basalganglia_/Brain_Putamen_basal_ganglia_/g
s/Brain_Spinalcord_cervicalc_1_/Brain_Spinal_cord_cervical_c_1_/g
s/Brain_Substantianigra/Brain_Substantia_nigra/g
s/Breast_MammaryTissue/Breast_Mammary_Tissue/g
s/Cells_EBV_transformedlymphocytes/Cells_EBV_transformed_lymphocytes/g
s/Cells_Transformedfibroblasts/Cells_Cultured_fibroblasts/g
s/Esophagus_GastroesophagealJunction/Esophagus_Gastroesophageal_Junction/g
s/FallopianTube/Fallopian_Tube/g
s/Heart_AtrialAppendage/Heart_Atrial_Appendage/g
s/Heart_LeftVentricle/Heart_Left_Ventricle/g
s/MinorSalivaryGland/Minor_Salivary_Gland/g
s/Skin_NotSunExposed_Suprapubic_/Skin_Not_Sun_Exposed_Suprapubic_/g
s/Skin_SunExposed_Lowerleg_/Skin_Sun_Exposed_Lower_leg_/g
s/SmallIntestine_TerminalIleum/Small_Intestine_Terminal_Ileum/g
s/WholeBlood/Whole_Blood/g
EOF

cat > ~/tmp/merge.cummings.w.gtex.tpm.sh <<'EOF'
#$ -S /bin/bash
#$ -o /net/home/nconnally/data/coding_variant_selection_cummings_2019/eo_files/
#$ -e /net/home/nconnally/data/coding_variant_selection_cummings_2019/eo_files/
#$ -V
#$ -l h_vmem=20G



cd /net/home/nconnally/data/coding_variant_selection_cummings_2019/selected_snps

awk '
BEGIN{print "SNP\tsymbol\tensg\tchr\tpos\ta1\ta2\ttissue\ttranscript_per\ttotal_tpm"}
NR == FNR {
tpm[$1,$3]=$4
}
NR != FNR {
print $0 "\t" tpm[$3,$8]
}' \
/net/home/nconnally/data/gtex_v8_tpm/data/gtex_v8.median_tpm.melted.tsv \
<(sed -rf cummings.gtex.tissue.key.sh ukbb.coding.snps.cummings.2019.all_tissues.0.1.tsv) \
> ukbb.coding.snps.cummings.w.gtex.tpm.tsv
EOF
