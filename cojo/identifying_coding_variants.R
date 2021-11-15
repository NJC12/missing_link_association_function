library(dplyr)
options(tibble.width = Inf)
library(ggplot2)
library(reshape2)

## Uses data downloaded in download_cummings_2019_data.sh

tab <- read.table('~/pleiotropy/data/gtex/cummings_2019_w_gtex/ukbb.coding.snps.cummings.w.gtex.tpm.tsv', header=T, as.is=T, sep='\t')

tran <- tab %>%
  group_by(ensg) %>%
  mutate(fraction_max_tpm=total_tpm/max(total_tpm, na.rm=T))

## HDL and LDL
chol.tis <- c('Liver', 'Whole_Blood', 'Adipose_Subcutaneous', 'Adipose_Visceral_Omentum_')
chol.snps <- tran[tran$fraction_max_tpm >= 0.5 &
                  tran$transcript_per >= 0.25 &
                  tran$tissue %in% chol.tis &
                  !is.na(tran$fraction_max_tpm), ]

write.table(unique(chol.snps$SNP), '~/tmp/coding.snps.for.cholesterol.tsv', row.names=F, col.names=F, quote=F)

## Height
height.snps <- tran[tran$fraction_max_tpm >= 0.5 &
                  tran$transcript_per >= 0.25 &
                  tran$tissue == 'Muscle_Skeletal' &
                  !is.na(tran$fraction_max_tpm), ]

write.table(unique(height.snps$SNP), '~/tmp/coding.snps.for.height.tsv', row.names=F, col.names=F, quote=F)

## Type II diabetes

system('rm ~/tmp/coding.snps.for.t2d.chr*')

t2d.snps <- tran[tran$fraction_max_tpm >= 0.5 &
                  tran$transcript_per >= 0.25 &
                  tran$tissue == 'Pancreas' &
                  !is.na(tran$fraction_max_tpm), ]

for (chr in 1:22) {

  t2d.chr <- t2d.snps[t2d.snps$chr == chr, ]
  if (nrow(t2d.chr) == 0) next
  t2d.list1 <- paste(t2d.chr$chr, t2d.chr$pos, t2d.chr$a1, t2d.chr$a2, sep=':')
  t2d.list2 <- paste(t2d.chr$chr, t2d.chr$pos, t2d.chr$a2, t2d.chr$a1, sep=':')

  write.table(unique(c(t2d.list1, t2d.list2)), paste0('~/tmp/coding.snps.for.t2d.chr', chr, '.tsv'), row.names=F, col.names=F, quote=F)

}

ibd.tis <- c('Small_Intestine_Terminal_Ileum', 'Colon_Transverse', 'Colon_Sigmoid')
ibd.snps <- tran[tran$fraction_max_tpm >= 0.5 &
                 tran$transcript_per >= 0.25 &
                 tran$tissue %in% ibd.tis &
                 !is.na(tran$fraction_max_tpm), ]

for (chr in 1:22) {

  ibd.chr <- ibd.snps[ibd.snps$chr == chr, ]
  if (nrow(ibd.chr) == 0) next
  ibd.list1 <- paste(ibd.chr$chr, ibd.chr$pos, ibd.chr$a1, ibd.chr$a2, sep=':')
  ibd.list2 <- paste(ibd.chr$chr, ibd.chr$pos, ibd.chr$a2, ibd.chr$a1, sep=':')

  write.table(unique(c(ibd.list1, ibd.list2)), paste0('~/tmp/coding.snps.for.ibd.chr', chr, '.tsv'), row.names=F, col.names=F, quote=F)

}


bc.snps <- tran[tran$fraction_max_tpm >= 0.5 &
                  tran$transcript_per >= 0.25 &
                  tran$tissue == 'Breast_Mammary_Tissue' &
                  !is.na(tran$fraction_max_tpm), ]

for (chr in 1:22) {

  bc.chr <- bc.snps[bc.snps$chr == chr, ]
  if (nrow(bc.chr) == 0) next
  bc.list1 <- paste(bc.chr$chr, bc.chr$pos, bc.chr$a1, bc.chr$a2, sep=':')
  bc.list2 <- paste(bc.chr$chr, bc.chr$pos, bc.chr$a2, bc.chr$a1, sep=':')

  write.table(unique(c(bc.list1, bc.list2)), paste0('~/tmp/coding.snps.for.bc.chr', chr, '.tsv'), row.names=F, col.names=F, quote=F)

}
