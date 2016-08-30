library(DESeq2)
library(sva)

perc_slight_diff = 1
min_strong_diff = 1 + as.double(commandArgs(T)[1])
set.seed(1463)

counts = read.table("dKD_counts", row.names=1)
colData = data.frame(contrast = c("c1","c2","c3"))
names(counts) = colData$contrast
dds = DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~contrast)
dds = DESeq(dds)
mean_disp_dKD=data.frame(mean=mcols(dds)$baseMean, disp=mcols(dds)$dispersion)
mean_disp_dKD$disp[is.na(mean_disp_dKD$disp)]=10
mean_disp_dKD$mean[mean_disp_dKD$mean==0]=0.1

counts = read.table("dKD_SMG7_counts", row.names=1)
colData = data.frame(contrast = c("c1","c2","c3"))
names(counts) = colData$contrast
dds = DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~contrast)
dds = DESeq(dds)
mean_disp_SMG7=data.frame(mean=mcols(dds)$baseMean, disp=mcols(dds)$dispersion)
mean_disp_SMG7$mean[is.na(mean_disp_SMG7$disp)]=0.1
mean_disp_SMG7$disp[is.na(mean_disp_SMG7$disp)]=10
mean_disp_SMG7$mean[mean_disp_SMG7$mean==0]=0.1
disp_fun = dispersionFunction(dds)
disp_coefs=as.double(attr(disp_fun,"coefficients"))

obs_log2FC = read.table("log2FCs_moderated",row.names=1)
obs_log2FC[is.na(obs_log2FC)] = 0

SMG6_fc = 2^(obs_log2FC[,1]*0.57)
SMG7_fc = 2^(obs_log2FC[,2])

# sligthly different genes
n_genes = length(mean_disp_SMG7$mean)
rand_index = which(runif(n_genes) < perc_slight_diff)
rand_sign = sample(x=c(1,-1),size=length(rand_index),replace=T)
rand_diff = rep(1,n_genes)
rand_diff[rand_index] = rand_diff[rand_index] * runif(length(rand_index),min=1,max=min_strong_diff)
rand_diff[rand_index][rand_sign==-1] = 1/rand_diff[rand_index][rand_sign==-1]
rand_diff = rand_diff * SMG6_fc

mean_disp = data.frame(genes=row.names(counts),mean=(mean_disp_SMG7$mean/SMG7_fc)*rand_diff)
mean_disp = mean_disp[order(mean_disp$mean,decreasing=T),]
mean_disp_SMG7 = mean_disp_SMG7[order(mean_disp_SMG7$mean,decreasing=T),]
mean_disp$disp = mean_disp_SMG7$disp
mean_disp = mean_disp[order(mean_disp$gene),]
mean_disp$genes = NULL

new_counts = read.table("scr_dKD_counts", row.names=1)
new_counts$new1=apply(mean_disp,1,function(x) rnbinom(n=1,mu=x[1],size=1/x[2])) # different definition of dispersion between deseq2 and rnbinom
new_counts$new2=apply(mean_disp,1,function(x) rnbinom(n=1,mu=x[1],size=1/x[2]))
new_counts$new3=apply(mean_disp,1,function(x) rnbinom(n=1,mu=x[1],size=1/x[2]))
colData = data.frame(legend = c("scr1","scr2","scr3","kd1","kd2","kd3","rescue1","rescue2","rescue3"))
names(new_counts) = colData$legend
trimLastChar = function (x) substr(x, 1, nchar(x) - 1)
colData$contrast = sapply(as.character(colData$legend), trimLastChar)
dds = DESeqDataSetFromMatrix(countData=new_counts, colData=colData, design=~contrast)

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

res_rescue = results(dds, contrast=c("contrast", "rescue", "kd"))
res_rescue$baseMean = NULL
res_rescue$stat = NULL
res_rescue$pvalue = NULL
res_rescue$lfcSE = NULL

write.table(res_rescue, file="results_simu", col.names=F, quote=F)
