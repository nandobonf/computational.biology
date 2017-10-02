#!/usr/bin/env Rscript

cat("### Loading required packages...\n")
# if (!require("optparse", quietly = T)) install.packages("optparse")
# library(optparse)
# suppressMessages(library(BiocInstaller))
# suppressMessages(if (!require('GWASTools', quietly = T)) biocLite('GWASTools'))
# suppressMessages(library(GWASTools))
# suppressMessages(if (!require('data.table', quietly = T)) install.packages('data.table'))
# suppressMessages(library(data.table))
# suppressMessages(library(tools, quietly = T))
# suppressMessages(library(stringr, quietly = T))

list.of.packages.cran <- c("tools", "stringr", "data.table", "optparse")
new.packages <- list.of.packages.cran[!(list.of.packages.cran %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages.cran, require, character.only = TRUE)
list.of.packages.bioc <- c("GWASTools")
new.packages <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
suppressMessages(source("https://bioconductor.org/biocLite.R"))
if(length(new.packages)) biocLite(new.packages)
cat("Loading required package: GWASTools\n")
suppressMessages(lapply(list.of.packages.bioc, require, character.only = TRUE))


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NA, 
              help="Text file with the columns P, CHROM and POS, it can also be compressed (.gz)", metavar="character"),
  make_option(c("-r", "--chromosome"), type="character", default="CHR", 
              help="Column name for the chromosome, [default = %default]", metavar="character"),
  make_option(c("-b", "--basepair"), type="character", default="BP", 
              help="Column name for the genomic position, [default = %default]", metavar="character"),
  make_option(c("-p", "--pvalue"), type="character", default="P", 
              help="Column name for the association p-value, [default = %default]", metavar="character"),
  make_option(c("-c", "--cut"), type="logical", default=F, 
              help="To speed up plotting, cut out snps with P > 0.05, [default = %default]", metavar="logical"),
  make_option(c("-m", "--manout"), type="character", default="out", 
              help="output file name [default = %default(.manhattan.png)]", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=5e-6, 
              help="threshold line for suggestive significance [default = %default]", metavar="numeric"),
  make_option(c("-q", "--qqout"), type="character", default="out", 
              help="output file name [default = %default(qq.png)]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


cat("### Reading input file...\n")
if (is.na(opt$file)) {stop("input files must be provided. See script usage (--help)")}
if (file_ext(opt$file) == 'gz') {
  data.man <- fread(paste0('gunzip -c ',opt$file), select = c(opt$pvalue, opt$chromosome, opt$basepair), data.table = F)
} else {
  data.man <- fread(paste(opt$file), select = c(opt$pvalue, opt$chromosome, opt$basepair), data.table = F)
}

cat("### Calculating genomic lambda...")
data.man <- data.man[!is.na(data.man[,1]),]
ntests = nrow(data.man)
lam <- round(x = median(qchisq(1-data.man[,1], 1))/qchisq(0.5, 1), 4)
cat("lambda =", paste(lam), "\n")
if (opt$cut == T) {
  data.man <- data.man[data.man[,1] <= 0.05,]
  cat("### Removing P values > 0.05 to increase plotting speed...\n")
}
data.man <- data.man[,!data.man[2] %in% c("X","x","Y","y","23","24")]
data.man[,c(2:3):= lapply(.SD, as.numeric), .SDcols = 2:3]
data.man <- data.man[order(data.man[,2], data.man[,3]),]
max <- ceiling(max(-log10(data.man[,1]))+1)
png(file = paste0(opt$manout, '.manhattan.png'), width = 3400, height = 2000, res = 300)
if (opt$cut == F) {
  manhattanPlot(data.man[,1], data.man[,2], ylim = c(0,max), signif = NULL)
} else {
  manhattanPlot(data.man[,1], data.man[,2], ylim = c(1,max), signif = NULL)
}
abline(h = -log10(opt$threshold), lty = 'dashed', col = 'blue')
abline(h = -log10(5e-8), lty = 'dashed', col = 'red')
invisible(dev.off())

cat("### Output files created:\n")
cat("### Manhattan plot written to:", paste0(opt$manout, '.manhattan.png\n'))

pvalues <- data.man[,1]
pvalues = sort.int(pvalues)
ypvs = -log10(pvalues)
xpvs = -log10(seq_along(ypvs)/ntests)
levels = as.integer((xpvs - xpvs[1])/(tail(xpvs, 1) - xpvs[1]) * 2000)
keep = c(TRUE, diff(levels) != 0)
levels = as.integer((ypvs - ypvs[1])/(tail(ypvs, 1) - ypvs[1]) * 2000)
keep = keep | c(TRUE, diff(levels) != 0)
keep = which(keep)
ypvs = ypvs[keep]
xpvs = xpvs[keep]
mx = head(xpvs, 1) * 1.05
my = max(mx * 1.15, head(ypvs, 1)) * 1.05
png(file = paste0(opt$qqout, '.qq.png'), res = 300, width = 1600, height = 1600)
plot(NA, NA, ylim = c(0, my), xlim = c(0, mx), xaxs = "i", 
     yaxs = "i", xlab = expression("- log"[10] * "(p-value), expected under null"), 
     ylab = expression("- log"[10] * "(p-value), observed"))
lines(c(0, mx), c(0, mx), col = "grey")
points(xpvs, ypvs, pch = 19, cex = 0.5)
title(main = paste0("Lambda = ", lam))
invisible(dev.off())
cat("### QQ-plot plot written to:", paste0(opt$qqout, '.qq.png\n'))
