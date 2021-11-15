Much of the COJO code is identifying coding SNPs and preparing the TOPMed data used as an LD reference.
- First run download_cummings_2019_data.sh, which downloads data on expression from Cummings et al. 2019
- Then run identifying_coding_variants.R, which filters out coding variants by how expressed they are in relevant tissues and how commonly they appear in exons of the different splice variants in that tissue (see the methods of the paper for more details).
- preparing_topmed_ld_reference.sh extracts 4,700 Europeans subjects from TOPMed to use as our reference LD
- [pheno]_topmed_reference.sh pulls out the SNPs used for each specific GWAS (this step includes a liftover to make builds match).
- run_cojo_[pheno].sh actually performs COJO (in addition to some data munging to make the relevant files, and some liftover).
