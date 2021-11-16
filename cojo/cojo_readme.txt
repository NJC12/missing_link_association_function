Much of the COJO code is identifying coding SNPs and preparing the TOPMed data used as an LD reference.
- First run download_cummings_2019_data.sh, which downloads data on expression from Cummings et al. 2019
- Then run identifying_coding_variants.R, which filters out coding variants by how expressed they are in relevant tissues and how commonly they appear in exons of the different splice variants in that tissue (see the methods of the paper for more details).
- preparing_topmed_ld_reference.sh extracts 4,700 Europeans subjects from TOPMed to use as our reference LD
- [pheno]_topmed_reference.sh pulls out the SNPs used for each specific GWAS (this step includes a liftover to make builds match).
- run_cojo_[pheno].sh actually performs COJO (in addition to some data munging to make the relevant files, and some liftover).

This code depends on the following datasets.
- GWAS summary statistics from
- Zhang et al. 2020 (http://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/gwas-summary-associations-breast-cancer-risk-2020/)
- Liu et al. 2015 and Goyette et al. 2015 (https://www.ibdgenetics.org/downloads.html or ftp://ftp.sanger.ac.uk/pub/consortia/ibdgenetics/iibdgc-trans-ancestry-filtered-summary-stats.tgz)
- Mahajan et al. 2018 (http://diagram-consortium.org/downloads.html)
  - Specifically, statistics "Adjusted for BMI" were used.
- Our conditioned analysis of the above datasets (https://doi.org/10.5061/dryad.612jm644q)
- Our GWAS for height, LDL, and HDL (https://doi.org/10.5061/dryad.612jm644q)
