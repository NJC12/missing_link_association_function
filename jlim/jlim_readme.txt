N.b. A newer version of JLIM is likely coming soon that would remove the need to use genotypic data, instead using summary statistics. This would allow JLIM to be run much more simply than in this script. If it is available, USE THE NEWER JLIM, as it will be much easier.

This code depends on several tools.
- bcftools: https://samtools.github.io/bcftools/bcftools.html
- vcftools: https://vcftools.github.io/index.html
- R: https://www.r-project.org/
- JLIM: https://github.com/cotsapaslab/jlim
- Perl: https://www.perl.org/

This code depends on the following datasets.
- GTEx: Summary statistics are freely available (https://gtexportal.org/home/datasets). Genotypes can only be used with permission of the consortium (https://gtexportal.org/home/protectedDataAccess), or once datasets are posted on dbGAP, with permission of the NIH (https://gtexportal.org/home/documentationPage).
- GWAS summary statistics from
- Zhang et al. 2020 (http://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-associations-breast-cancer-risk-2020/)
- Liu et al. 2015 and Goyette et al. 2015 (https://www.ibdgenetics.org/downloads.html or ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz)
- Mahajan et al. 2018 (http://diagram-consortium.org/downloads.html)
  - Specifically, statistics "Adjusted for BMI" were used.
- Our GWAS for height, LDL, and HDL (https://doi.org/10.5061/dryad.612jm644q)

Variables
- Several global variables must be set
- Window size: The distance around each GWAS peak to check for the TSS.
- Genotype: The location of the GTEx genotype data.
- Several variables are passed (see examples below)
- dir: The parent directory for the analysis (structure shown below)
- prefix: The string added to the beginning of many generated file names in order to desparately preserve some semblance of order.
- gwas_sum_stats: The file containing your GWAS summary statistics.
- tissue: The name of the GTEx tissue you are analyzing.
- ethnicIDs: With current GTEx sample sizes, this analysis is likely to be severely underpowered on non-European ethnic groups. This text file lists the IDs of the GTEx samples you would like to include in the analysis. Despite the name, this sample can also be used to perform a sex-specific GWAS (e.g. on breast cancer).
- covar: Covariates for the GTEx subjects included.

Example passed variables
| directory  | /home/nconnally/pleiotropy_equilibrium/jlim/ldl                                                                                   |
| prefix     | ukbb_ldl                                                                                                                          |
| gwas_stats | /home/nconnally/pleiotropy_equilibrium/jlim/ldl_quantile/raw_data/ukbb.ldl.quantile.all.chr.output.conditioned.cma.cojo           |
| tissue     | liver                                                                                                                             |
| ethnic_ids | /net/equilibrium/dbGaP/analysis/phs000424/GTEx.v7/Liver/Liver-Covariates/Whites/Liver_Europeans.txt                               |
| expr       | /net/equilibrium/dbGaP/analysis/phs000424/GTEx.v7/GTEx_Analysis_v7_eQTL_expression_matrices/Liver.v7.normalized_expression.bed.gz |
| covariates | /net/equilibrium/dbGaP/analysis/phs000424/GTEx.v7/GTEx_Analysis_v7_eQTL_covariates/Liver.v7.covariates.txt                        |

Sections
- By flipping the variables in "Which sections to run" between 0 and 1, different sections can be run independently.

Directory structure (within the parent directory defined by "$dir")
- tmp_src: Where the main program (run_jlim.sh) is stored. Simply referred to as temporary because I was actually running this whole operation out of an emacs org-mode file.
- eo_files: For error and output files generated in running the scripts created by this program.
- main_pheno_peaks: When the R program 2019.10.03.call.peaks.R is run, it deposits the blocks of summary statistics around each genome-wide-significant peak from the GWAS here.
- main_pheno_peaks_no_dups: As above, but removes duplicated SNPs. The fact that this process is separate from the above is a superfluity, but there isn't enough selection pressure to remove it.
- secondary_geno_peaks: When the relevant regions are selected from the GWAS, the GTEx genoytpes for the same regions are extracted and stored here.
- tmp_qsub_files: This program submits separate jobs for different chromosomes as a form of NSPP (Noah's stupid parallel processing). Files for these different jobs are stored in this directory so they can be easily deleted later.
- tmp_vcf_files: Regions of the genome extracted from GWAS are stored as VCFs here, but are quickly altered and moved so as not to take up so much hard drive space.
- secondary_phenotypes_$tissue: The secondary phenotype (i.e. gene expression for this particular tissue) is stored here. Because one GWAS phenotype will often be compared to expression in multiple tissues, several such directories will often coexist.
- covariates: This covariates are taken from the GTEx data, manipulated into a more useful form, and stored here.
- secondary_assoc_perm: eQTL calling is performed using the GTEx gene expression data and GTEx genotypes. The results are stored here.
- reference_ld: JLIM generates reference LD files, which are stored here.
- jlim_cfg_$tissue: Stores intermediary results for JLIM.
- jlim_out_$tissue: Stores the actual JLIM results.
