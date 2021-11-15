The Coloc package implements the method first described in Giambartolomei et al. 2013.
- https://cran.r-project.org/web/packages/coloc/

For ineffeciency's sake, I wrote a separate small script for analyzing each phenotype.
- Each script extract_[pheno]_gwas_regions.sh breaks the GWAS summary stats into separate files based on the indices generated in call_peaks.
- Each script extract_[pheno]_gtex_[tissue]_regions.sh does the same for the GTEx summary statistics.
- Then, run_calc_coloc_[pheno]_[tissue].sh generates a file to submit each job to calc_coloc.R and then runs the file.

For the GWAS which was ran from UKBB data, there is also a script to extract the minor allele frequencies (MAFs), which are needed for coloc.
- calc_ukbb_mafs.sh
