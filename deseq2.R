library(DESeq2)
library(ggplot2)
library(dplyr)
library(sva)

counts = read.table("counts", row.names=1)
colData = data.frame(legend = c("scr1","scr2","scr3","kd1","kd2","kd3","rescue1","rescue2","rescue3"))
trimLastChar = function (x) substr(x, 1, nchar(x) - 1)
colData$contrast = sapply(as.character(colData$legend), trimLastChar)
names(counts) = colData$contrast
dds = DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~contrast)

# sva
dds = estimateSizeFactors(dds)
dat = counts(dds, normalized=TRUE)
dat = dat[rowMeans(dat)>5,]
mod = model.matrix(~contrast,colData(dds))
mod0 = model.matrix(~1,colData(dds))
svseq = svaseq(dat, mod, mod0, n.sv=1)
dds$sv = svseq$sv
design(dds) = ~ sv + contrast

dds = DESeq(dds)

# deseq
normCounts = counts(dds, normalized=T)
colnames(normCounts) = colData$legend
res_kd = results(dds, contrast=c("contrast", "kd", "scr"))
res_kd$gene_id = row.names(res_kd)
res_rescue = results(dds, contrast=c("contrast", "rescue", "kd"))
res_rescue$gene_id = row.names(res_rescue)

# output
res = merge(as.data.frame(res_kd), as.data.frame(res_rescue), by="gene_id")
write.table(res, file="results", col.names=F, row.names=F, quote=F)
