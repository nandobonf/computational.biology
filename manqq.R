#!/usr/bin/env Rscript


if (!require("optparse", quietly = T)) install.packages("optparse")
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NA, 
              help="Text file with the columns P, CHROM and POS, it can also be compressed (.gz)", metavar="character"),
  make_option(c("-", "--chromosome"), type="character", default="CHR", 
              help="Column name for the chromosome, [default = %default]", metavar="character"),
  make_option(c("-b", "--basepair"), type="character", default="BP", 
              help="Column name for the genomic position, [default = %default]", metavar="character"),
  make_option(c("-p", "--pvalue"), type="character", default="P", 
              help="Column name for the association p-value, [default = %default]", metavar="character"),
  make_option(c("-c", "--cut"), type="logical", default=F, 
              help="To speed up plotting, cut out snps with P > 0.05, [default = %default]", metavar="logical"),
  make_option(c("-m", "--manout"), type="character", default="out", 
              help="output file name [default = %default(.manhattan.png)]", metavar="character"),
  make_option(c("-q", "--qqout"), type="character", default="out", 
              help="output file name [default = %default(qq.png)]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#inputname <- "~/UKB/icd10/final.results/UKB.results.IBS.icd10.logistic.annotated.gz" # only for testing
suppressMessages(library(BiocInstaller))
suppressMessages(if (!require('GWASTools', quietly = T)) biocLite('GWASTools'))
suppressMessages(library(GWASTools))
suppressMessages(if (!require('data.table', quietly = T)) install.packages('data.table'))
suppressMessages(library(data.table))
suppressMessages(library(tools, quietly = T))
suppressMessages(library(stringr, quietly = T))

# outfilename <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(inputname))

if (is.na(opt$file)) {stop("input files must be provided. See script usage (--help)")}
if (file_ext(opt$file) == 'gz') {
  data.man <- fread(paste0('gunzip -c ',opt$file), select = c(opt$pvalue, opt$chromosome, opt$basepair), data.table = F)
} else {
  data.man <- fread(paste(opt$file), select = c(opt$pvalue, opt$chromosome, opt$basepair), data.table = F)
}

if (opt$cut == T) {
  data.man <- data.man[data.man[,1] <= 0.05,]
}
data.man <- data.man[!is.na(data.man[,1]),]
data.man <- data.man[order(data.man[,2], data.man[,3]),]
max <- ceiling(max(-log10(data.man[,1]))+1)
png(file = paste0(opt$manout, '.manhattan.png'), width = 3400, height = 2000, res = 300)
manhattanPlot(data.man[,1], data.man[,2], ylim = c(0,max), signif = NULL)
abline(h = -log10(5e-6), lty = 'dashed', col = 'blue')
abline(h = -log10(5e-8), lty = 'dashed', col = 'red')
dev.off()
data.man <- as.data.frame(data.man)
lam <- round(x = median(qchisq(1-data.man[,1][!is.na(data.man[,1])],1))/qchisq(0.5,1), 4)
png(file = paste0(opt$qqout, '.qq.png'), res = 300, width = 1600, height = 1600)
qqPlot(data.man[,1], ci = F, main = paste('lambda = ', lam))
dev.off()
