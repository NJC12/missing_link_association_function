#!/bin/bash

dir="${vars[0]}"
prefix="${vars[1]}"
gwas_sum_stats="${vars[2]}"
tissue="${vars[3]}"
ethnicIDs="${vars[4]}"
expr="${vars[5]}"
covar="${vars[6]}"

mkdir -p $dir/eo_files $dir/tmp_src

# This was added because of problems with the cluster on which the code was being run.
# It should not be necessary for others to include.
TERM=xterm-color
export TERM

#######################################################
#                Set global variables                 #
#######################################################

# The size on either side of the lead SNP to extract. Final size will be windowsize*2+1
windowsize=100000

# The vcf with the genotypes for GTEx.v7
genotype="/net/equilibrium/dbGaP/data_downloads/projects/proj13919/57936/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz"

# List of transcription start sites for each gene
tss="/net/equilibrium/dbGaP/analysis/phs000424/GTEx.v7/gene_annotations/bioMart_genes.hg19.bed"

# The range around each window to check for TSSs
tss_range=1000000

#######################################################
#                Which sections to run                #
#######################################################

# Set to 0 to skip, 1 to run.
# Generally, once a section is run all the following sections should be too.

call_peaks=1
# Subsection of call_peaks
no_x=1
uniq_tr1_pos=1
# Subsection of extract_from_vcf
# Should not normally need to be turned on
index_vcf=0
extract_from_vcf=1
find_proximal_genes=1
fix_covar=1
assoc_perm=1
create_ld=1
jlim_gencfg=1
jlim=1
collect_jlim=1


#######################################################
# Gets the R module working after environment changed #
#######################################################

MODULEPATH=$MODULEPATH:/net/module/Modules/versions:/net/module/Modules/3.2.10/modulefiles:/net/module/Modules/modulefiles
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/net/module/sw/tcl/8.5.19/lib
export MODULEPATH
export LD_LIBRARY_PATH

#######################################################
#         Calling peaks in the main phenotype         #
#######################################################

if [ $call_peaks -gt 0 ]; then

  mkdir -p ${dir}/main_pheno_peaks

  # The arguments are
  # Input summary stats
  # Output files base name
  # Window size on either side
  # Column with chromosome
  # Column with position
  # Column with variant ID
  # Column with p-value
  # Column with number of samples
  # Column with minor allele frequency
  module load R
  Rscript /home/nconnally/pleiotropy/src/jlim/code/2019.01.03.call.peaks.R \\
          $gwas_sum_stats \\
          ${dir}/main_pheno_peaks/${prefix} \\
          $windowsize \\
          1 \\
          3 \\
          2 \\
          13 \\
          9 \\
          5

fi

if [ $no_x -gt 0 ]; then

  if [ -f "${dir/main_pheno_peaks/WITH_X.${prefix}.indexSNP.tsv" ]; then
      awk '$1 != "X"' ${dir}/main_pheno_peaks/WITH_X.${prefix}.indexSNP.tsv > ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv
  else
      mv ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv ${dir}/main_pheno_peaks/WITH_X.${prefix}.indexSNP.tsv
      awk '$1 != "X"' ${dir}/main_pheno_peaks/WITH_X.${prefix}.indexSNP.tsv > ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv
  fi
fi


#######################################################
#         Calling peaks in the main phenotype         #
#######################################################

if [ $uniq_tr1_pos -gt 0 ]; then

  mkdir -p ${dir}/main_pheno_peaks_no_dups/

  while read -r line; do

      chr=$(echo $line | awk '{print $1}')
      start=$(echo $line | awk '{print $4}')
      end=$(echo $line | awk '{print $5}')

      awk '{if (pos[$2] != 1) {print $0; pos[$2] = 1}}' ${dir}/main_pheno_peaks/${prefix}.${chr}.${start}.${end}.txt \
          > ${dir}/main_pheno_peaks_no_dups/${prefix}.${chr}.${start}.${end}.txt

  done < ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv
fi


#######################################################
#        Extract the same peaks from the vcf          #
#######################################################

if [ $index_vcf -gt 0 ]; then

  mkdir -p ${dir}/secondary_geno_peaks

  # awk 'NR > 1 {print $1"\t"$4"\t"$5}' ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv > \
  #     ${dir}/secondary_geno_peaks/bcf.regions.tsv

  module load bcftools

  bcftools index -c -f $genotype

fi

if [ $extract_from_vcf -gt 0 ]; then

  mkdir -p ${dir}/secondary_geno_peaks

  # We have all the peaks we need stored in the indexSNP file
  # but we need to convert it into a format that can be used to extract regions from the vcf.
  # vcftools .bed files also need a header

  mkdir -p ${dir}/tmp_qsub_files
  mkdir -p ${dir}/tmp_vcf_files

  # This builds up a string of all the jobs we have to wait for the completion of
  job_ids=''

  while read -r line; do
      chr=$(echo $line | awk '{print $1}')
      start=$(echo $line | awk '{print $4}')
      end=$(echo $line | awk '{print $5}')

      # Skip the header row
      if [ "$chr" == "CHR" ]; then
          continue
      fi

      DATE=`date '+%Y-%m-%d %H:%M'`

      cat > ${dir}/tmp_qsub_files/bcftools.${chr}.${start}.${end}.sh <<EOF2
#!/bin/bash
#$ -e /home/nconnally/pleiotropy/src/jlim/experiments/eo/bcftools.err
#$ -o /home/nconnally/pleiotropy/src/jlim/experiments/eo/bcftools.out

# Generated by \$0 on \$DATE

# Gets the modules working after environment changed

MODULEPATH=\$MODULEPATH:/net/module/Modules/versions:/net/module/Modules/3.2.10/modulefiles:/net/module/Modules/modulefiles
LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/net/module/sw/tcl/8.5.19/lib
export MODULEPATH
export LD_LIBRARY_PATH

# Actually does the work I want done

module load bcftools
module load vcftools

# My list of European samples included several that were removed in qc, so I use
# --force-samples to make bcftools ignore them.

bcftools view $genotype \\
       --regions ${chr}:${start}-${end} \\
       --samples-file $ethnicIDs \\
       --force-samples \\
       --output-type z \\
       --output-file ${dir}/tmp_vcf_files/${prefix}.${chr}.${start}.${end}.vcf.gz

vcftools --gzvcf ${dir}/tmp_vcf_files/${prefix}.${chr}.${start}.${end}.vcf.gz \
       --plink \
       --out ${dir}/secondary_geno_peaks/${prefix}.${chr}.${start}.${end}       

EOF2

      # Each time a job is submitted, its ID is added to this list of job IDs.
      # This is so that we can run a dummy job below to wait for all of these to finish.
      job_ids=${job_ids},$(qsub -V -terse ${dir}/tmp_qsub_files/bcftools.${chr}.${start}.${end}.sh)

  done < ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv

  # This job is just to make qsub wait for the other jobs to finish before proceding.
  # -hold_jid waits for the list of IDs that we have built by recording each job we submit.
  # -sync y prevents the script from moving on until this job is finished.
  # -b y tells us that our job is a command, not a script. In this case it's echo.
  qsub -hold_jid ${job_ids} -sync y -b y echo 'Done'

fi


#######################################################
#         Find genes proximal to pheno peaks          #
#######################################################

if [ $find_proximal_genes -gt 0 ]; then

  module load R

  mkdir -p ${dir}/secondary_phenotypes_$tissue

  # Transcription start sites (TSS) are stored in
  # /net/equilibrium/dbGaP/analysis/phs000424/GTEx.v7/gene_annotations/bioMart_genes.hg19.bed

  # Arguments are:
  # Peaks (the indexSNP file)
  # Range around each window to select genes from
  # Output file
  # Expression data that will be written into a phenotype file
  Rscript /home/nconnally/pleiotropy/src/jlim/code/2019.01.10.extract.genes.near.peaks.R \\
          ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv \\
          $tss_range \\
          ${dir}/secondary_phenotypes_${tissue}/${prefix}.${tissue} \\
          $expr

fi


#######################################################
#              Manipulated covariate file             #
#######################################################

if [ $fix_covar -gt 0 ]; then

  mkdir -p ${dir}/covariates

  # Transpose covariate
  # And replicate the first column so that we have FID and IID
  # This should probably be a script of its own, but it's just a long awk command that
  # enters everything into a matrix, then prints it out transposed.
  awk '{ for (i=1; i<=NF; i++){
           cells[NR][i]=$i }} END {
           for(i=1; i<=NF; i++){
             printf cells[1][i]"\t"; for(j=1; j<NR; j++){
           printf cells[j][i]"\t" }; print cells[NR][i] }}' $covar | \\
               sed '1s/ID ID/FID IID/' \\
                   > ${dir}/covariates/${prefix}.${tissue}.covariates.txt

fi


#######################################################
#           Linear regression & permutation           #
#######################################################

if [ $assoc_perm -gt 0 ]; then

  mkdir -p ${dir}/secondary_assoc_perm_$tissue
  mkdir -p ${dir}/secondary_geno_peaks/dups_to_remove
  mkdir -p ${dir}/tmp_qsub_files/

  # This builds up a string of all the jobs we have to wait for the completion of
  job_ids=''

  while read -r line; do
      chr=$(echo $line | awk '{print $1}')
      start=$(echo $line | awk '{print $4}')
      end=$(echo $line | awk '{print $5}')

      # Skip the header row and X chromosome
      if [ "$chr" == "CHR" ]; then
          continue
      fi
      if [ "$chr" == "X" ]; then
          continue
      fi

      plink_file="${dir}/secondary_geno_peaks/${prefix}.${chr}.${start}.${end}"
      pheno_file="${dir}/secondary_phenotypes_${tissue}/${prefix}.${tissue}.phenos.${chr}.${start}.${end}.txt"
      n_cols=$(awk 'NR == 1 {print NF}' $pheno_file)

      if [ -f plink_file ]; then
          echo "THE FILE EXISTS: ${plink_file}"
      fi
      if [ "$chr" == "X" ]; then
          echo "Ignoring X chromosome"
      fi

      # There are multiple variants at the same position in many of these files, which causes plink to throw an error.
      # I remove all duplicates, leaving no copies of each.
      # This lists the variants that are removed with --exclude below.

      awk '{if(seen[$4] == 1) {print $1" "$4" "$4+1" tmp"} else {seen[$4] = 1}}' ${plink_file}.map > \\
          ${dir}/secondary_geno_peaks/dups_to_remove/${prefix}.${chr}.${start}.${end}.dups.to.remove.range.txt

      new_plink_file="${dir}/secondary_geno_peaks/${prefix}.no.dups.${chr}.${start}.${end}"

      # JLIM requires all the assoc and perm files to be in a directory of their own.
      # So we create the directory and move them there.
      mkdir -p ${dir}/secondary_assoc_perm_${tissue}/locus.${chr}.${start}.${end}

      cat > ${dir}/tmp_qsub_files/plink.${chr}.${start}.${end}.sh <<EOF2
#!/bin/bash
#$ -e /home/nconnally/pleiotropy/src/jlim/experiments/eo/plink.err
#$ -o /home/nconnally/pleiotropy/src/jlim/experiments/eo/plink.out
#$ -V

# I have to use a MAF filter, because otherwise JLIM filters them out of one file, but not the other, and throws an error
# If I only use the built-in plink functions to remove triallelic SNPs, it won't find ones where the different options have different SNP names
# If --missing-code is not set, plink interprets only NA as missing, while our data have 0's.
/home/nconnally/bin/plink --file $plink_file \\
                        --exclude range ${dir}/secondary_geno_peaks/dups_to_remove/${prefix}.${chr}.${start}.${end}.dups.to.remove.range.txt \\
                        --maf 0.05 \\
                        --snps-only \\
                        --biallelic-only strict \\
                        --geno 0.00 \\
                        --missing-code -9,0,NA,na \\
                        --recode \\
                        --out $new_plink_file

# Do use --snps-only and --biallelic-only strict in order to only have one variant at the same site
# (latter may be unnecessary, but doesn't hurt).
# Using hide-covar with --linear keeps the regression coefficient for each SNP from being printed
# with each variant.
# --freq produces the MAFs I need for coloc
/home/nconnally/bin/plink --file $new_plink_file \\
                        --pheno $pheno_file \\
                        --all-pheno \\
                        --linear mperm=1000 hide-covar \\
                        --covar ${dir}/covariates/${prefix}.${tissue}.covariates.txt \\
                        --no-const-covar \\
                        --seed 19930224 \\
                        --mperm-save-all \\
                        --allow-no-sex \\
                        --freq \\
                        --out ${dir}/secondary_assoc_perm_${tissue}/locus.${chr}.${start}.${end}/${prefix}.${tissue}.${chr}.${start}.${end}

if [ ! -f ${dir}/secondary_assoc_perm_${tissue}/locus.${chr}.${start}.${end}/ldsc_hub.ped ]; then
  cp ${dir}/secondary_geno_peaks/${prefix}.no.dups.${chr}.${start}.${end}.ped ${dir}/secondary_assoc_perm_${tissue}/locus.${chr}.${start}.${end}/${prefix}.ped
fi

gzip -fq ${dir}/secondary_assoc_perm_${tissue}/locus.${chr}.${start}.${end}/*

EOF2

      # Each time a job is submitted, its ID is added to this list of job IDs.
      # This is so that we can run a dummy job below to wait for all of these to finish.
      job_ids=${job_ids},$(qsub -V -terse ${dir}/tmp_qsub_files/plink.${chr}.${start}.${end}.sh)

  done < ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv

  # This job is just to make qsub wait for the other jobs to finish before proceding
  # This job is just to make qsub wait for the other jobs to finish before proceding.
  # -hold_jid waits for the list of IDs that we have built by recording each job we submit.
  # -sync y prevents the script from moving on until this job is finished.
  # -b y tells us that our job is a command, not a script. In this case it's echo.
  qsub -hold_jid ${job_ids} -sync y -b y echo 'Done'

fi


#######################################################
#                   Create LD files                   #
#######################################################

if [ $create_ld -gt 0 ]; then

  mkdir -p ${dir}/reference_ld

  module load perl
  module load vcftools

  # Sung has perl code for pulling out the LD.
  # It does not run very quickly.

  # The arguments are
  # 1KG data
  # My index SNP file, which I have used in many of the previous sections
  # The output folder

  # The index file HAS TO BE A TSV. Space-separated does not work.
  # Need to add the directory with tabix to my path
  PATH=$PATH:/net/equilibrium/dbGaP/imputation/bin/

  perl /home/nconnally/pleiotropy/src/jlim/code/2018.10.24.noah.edit.fetch.refld0.EUR.pl \\
       /net/equilibrium/dbGaP/post-imputation/1KG \\
       ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv \\
       ${dir}/reference_ld

fi


#######################################################
#                       JLIM CFG                      #
#######################################################

if [ $jlim_gencfg -gt 0 ]; then

  module load R
  # # Also have to install the JLIM package
  # Rscript -e 'install.packages("getopt", "~/bin/R", repos="http://cran.r-project.org")' 
  # R CMD INSTALL -l ~/bin/R /home/nconnally/pleiotropy/src/jlim/code/jlim/jlimR_1.0.2.tar.gz
  # And tell R to read in packages from this directory
  export R_LIBS=/home/nconnally/bin/R:$R_LIBS

  # And a new directory for each output
  mkdir -p ${dir}/jlim_cfg_$tissue

  bash /home/nconnally/pleiotropy/src/jlim/code/jlim/bin/jlim_gencfg.sh \\
       --tr1-name ${prefix} \\
       --tr1-dir ${dir}/main_pheno_peaks_no_dups \\
       --tr2-dir ${dir}/secondary_assoc_perm_${tissue} \\
       --idxSNP-file ${dir}/main_pheno_peaks/${prefix}.indexSNP.tsv \\
       --refld-dir ${dir}/reference_ld \\
       --out ${dir}/jlim_cfg_${tissue}/${prefix}.cfg.tsv

fi


#######################################################
#                         JLIM                        #
#######################################################

if [ $jlim -gt 0 ]; then

  module load R
  export R_LIBS=/home/nconnally/bin/R:$R_LIBS

  mkdir -p ${dir}/jlim_out_$tissue

  bash /home/nconnally/pleiotropy/src/jlim/code/jlim/bin/run_jlim.sh \\
       ${dir}/jlim_cfg_${tissue}/${prefix}.cfg.tsv \\
       0.8 \\
       ${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.out.tsv \\
       > ${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.out.log

fi


#######################################################
#                    Collect JLIM                     #
#######################################################

if [ $collect_jlim -gt 0 ]; then

  # This command writes the cfj and jlim output files together
  awk '{if (NR==FNR){snp[NR]=$3; chr[NR]=$2; pos[NR]=$4; split($14,gene1,"ENSG"); split(gene1[2],gene2,".assoc"); gene[NR]="ENSG"gene2[1]} else{stats[FNR]=$0}} END{for(i=1; i<=FNR; i++) print gene[i]"\t"snp[i]"\t"chr[i]"\t"pos[i]"\t"stats[i]}' \\
${dir}/jlim_cfg_${tissue}/${prefix}.cfg.tsv \\
${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.out.tsv > \\
${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.table.tsv

  module load R
  # This Rscript adds the gene names
  # The arguments are:
  # The JLIM table generated by the awk command above
  # A file with both ensembl gene and transcript names (you shouldn't have to change this)
  # A file with both ensembl gene names and regular gene names (you shouldn't have to change this)
  # The file to write the output to (here I'm just writing to the same one.
  Rscript /home/nconnally/pleiotropy/src/jlim/code/2019.03.22.add.gene.names.to.table.R \\
          ${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.table.tsv \\
          /home/nconnally/pleiotropy/data/reference/ensGene.txt.gz \\
          /home/nconnally/pleiotropy/data/reference/ensemblToGeneName.txt.gz \\
          ${dir}/jlim_out_${tissue}/${prefix}.${tissue}.jlim.table.tsv

fi
