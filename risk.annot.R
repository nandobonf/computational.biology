## define regions in META SELF 
library(plyr)
P_threshold <- 5e-6
d <- subset(repli, `P-value` <= P_threshold)
d <- data.table(chrom = d$CHR
                , name = as.character(d$SNP)
                , start = d$BP-250000
                , stop = d$BP+250000)
setorder(d, chrom, start)
d$start <- ifelse(d$start < 0, 0, d$start)
d <- ddply( d, "chrom", function(u) { 
  x <- c(NA, u$stop[-nrow(u)])
  y <- ifelse( is.na(x), 0, x )
  y <- cummax(y)
  y[ is.na(x) ] <- NA
  u$previous_stop <- y
  u
} )
d$new_group <- is.na(d$previous_stop) | d$start >= d$previous_stop
d$group <- cumsum( d$new_group )
d.final <-ddply( 
  d, .(chrom,group), summarize, 
  start=min(start), stop=max(stop), name=paste(name,collapse=", ")
)
d.final <- as.data.table(d.final)
d.final[,length_kb:=round(((stop-start)/1000),digits = 2)]

# add tag snp, p, or and se
d <- merge(d, repli, by.x = 2, by.y = 17, sort = F)
d <- as.data.table(d)
tags <- d[, .SD[which.min(`P-value`)], by = .(group)]
tags <- tags[, .(name, group, Effect, StdErr, Direction, `P-value`, Freq1, Allele1, Allele2)]
d.final <- merge(d.final, tags, by = 'group', sort = F)
colnames(d.final)<- c("locus_number","CHR","BP_start","BP_end","significant_SNPs","length_kb","TAG_SNP","BETA","SE","Direction","P","EAF", "EA", "OA")
d.final[,region_id:=paste0("chr",CHR,":",BP_start,"-", BP_end)]

# add mapping genes with XGR
library(XGR)
res <- xGR2nGenes(data = d.final$region_id, 
                  , format = "chr:start-end"
                  , decay.kernel = "rapid"
                  , distance.max = 0)
res <- as.data.table(res)
res.sum <- res[, .(mapping_genes = paste(Gene, collapse = ', ')),by = 'GR']
d.final <- merge(d.final, res.sum, by.x = "region_id", by.y = "GR", all.x = T)

# add eqtl genes GTEx v6
gtex <- readRDS("../Gtex_v6/significant_gtex_all.RDS")
gtex <- as.data.table(gtex)
head(gtex)
unique(gtex$gene_type)
gtex.cod <- gtex[gene_type %in% c("protein_coding","miRNA")]

d.final[,eqtl_genes:=NA]
for(i in 1:nrow(d.final)) {
  reg1 <- gtex.cod[snp_chrom == d.final$CHR[i] & snp_pos >= d.final$BP_start[i] & snp_pos <= d.final$BP_end[i], .(gene_name)]
  reg1 <- reg1[!duplicated(gene_name),]
  if(nrow(reg1) == 0){
    d.final$eqtl_genes[i] <- NA
  } else {
    d.final$eqtl_genes[i] <- paste(reg1$gene_name, collapse = ", ")
  }
}

all.genes <- unique(c(unlist(strsplit(d.final$mapping_genes, split = ", ")),unlist(strsplit(d.final$eqtl_genes, split = ", "))))
all.genes <- all.genes[!is.na(all.genes)]
fwrite(as.data.table(all.genes), "all.genes.IBS.meta.self.forGSEA", sep = "\t", col.names = F)

# MAP GENES WITH THE UCSC GENES
library(dplyr)
library(GenomicRanges)
library(Homo.sapiens)

d.final[,mapping_genes_UCSC:=NA]

for(i in 1:nrow(d.final)){
  mycoords.gr = d.final[i,.(CHR, BP_start, BP_end, region_id)] %>%
    mutate(CHR=paste0('chr', CHR)) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  print(mycoords.gr)
  
  maps <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr) %>% 
    as.data.frame %>% 
    merge(y = as.data.frame(org.Hs.egSYMBOL), by = "gene_id") %>%
    select_("symbol") %>%
    unlist(use.names = F) %>%
    paste0(collapse = " ") %>%
    gsub(pattern = " ", replacement = ", ")
  
  if(maps != "") {
    d.final$mapping_genes_UCSC[i] <- paste(maps)
  } else {
    d.final$mapping_genes_UCSC[i] <- NA}
}
