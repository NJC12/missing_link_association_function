#!/bin/bash
library(coloc)

gwas <- read.table(args[1], header=T, as.is=T)
gtex <- read.table(args[2], header=T, as.is=T)
outfile  <- args[3]
tis <- args[4]
n1 <- as.numeric(args[5])

gtex$chr <- sapply(gtex$variant_id, function(x) as.numeric(strsplit(x, '_')[[1]][1]))
gtex$bp <- sapply(gtex$variant_id, function(x) as.numeric(strsplit(x, '_')[[1]][2]))
gtex$maf <- as.numeric(gtex$maf)
gtex$slope <- as.numeric(gtex$slope)
gtex$slope_se <- as.numeric(gtex$slope_se)
gwas$n[gwas$n == 'nan'] <- NA
gwas$Chr <- as.numeric(gwas$Chr)
gwas$bp <- as.numeric(gwas$bp)
gwas$bC <- as.numeric(gwas$bC)
gwas$bC_se <- as.numeric(gwas$bC_se)
gwas$n <- as.numeric(gwas$n)
gwas$freq <- as.numeric(gwas$freq)

gtex.n <- mean(gtex$ma_count/(gtex$maf+0.0000001))
n2 <- round(mean(gtex.n[gtex.n > 0]))

calc.coloc <- function(sub.gwas, sub.gtex, gene) {
sub2.gtex <- sub.gtex[sub.gtex$gene_id == gene & sub.gtex$bp %in% sub.gwas$bp, ]
sub2.gwas <- sub.gwas[sub.gwas$bp %in% sub2.gtex$bp, ]
sub2.gwas$bC_se[sub2.gwas$bC_se <= 0] <- 0.01

coloc <- coloc.abf(dataset1=list(pvalues=sub2.gwas$pC,
                                 N=n1, MAF=sub2.gwas$freq, type="quant"),
                dataset2=list(pvalues=sub2.gtex$pval_nominal,
                              N=n2, MAF=sub2.gtex$maf, type="quant"))

return(coloc$summary)
}

call.coloc <- function(sub.gwas, sub.gtex, results) {
sub.gwas <- sub.gwas[!is.na(sub.gwas$bC), ]
sub.gtex <- sub.gtex[!is.na(sub.gtex$slope), ]

genes <- unique(sub.gtex$gene_id)

for (gene in genes) {
    coloc.summary <- calc.coloc(sub.gwas, sub.gtex, gene)
    h0 <- coloc.summary['PP.H0.abf']
    h1 <- coloc.summary['PP.H1.abf']
    h2 <- coloc.summary['PP.H2.abf']
    h3 <- coloc.summary['PP.H3.abf']
    h4 <- coloc.summary['PP.H4.abf']
    if (!is.na(h4) && h4 > results$h4[results$gene == gene]) {
       results$h0[results$gene == gene] <- h0
       results$h1[results$gene == gene] <- h1
       results$h2[results$gene == gene] <- h2
       results$h3[results$gene == gene] <- h3
       results$h4[results$gene == gene] <- h4
    }
}
return(results)
}

out <- data.frame(gene=unique(gtex$gene_id), h0=0, h1=0, h2=0, h3=0, h4=0, tissue=tis)

out <- call.coloc(gwas, gtex, out)

write.table(out, outfile, row.names=F, quote=F, sep='\t')
