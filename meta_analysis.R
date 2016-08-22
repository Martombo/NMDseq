# meta analysis functions
Fisher.test <- function(p) {
	Xsq <- -2*sum(log(p))
	p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
	return(p.val)
}

psumunif = function(x,n) {
	    fun = function(k) (-1)^k * choose(n,k) * ifelse(x > k,x-k,0)^(n)
		return(1/factorial(n) * sum(sapply(0:n, fun)))
}

tab=read.table("merge",header=T)
minFC=0.1
kd_cond=c("upf1_kd","smg6_kd","smg7_kd","dKD")
rescue_cond=c("upf1_rescue","smg6_rescue","smg7_rescue","dKD_smg6rescue","dKD_smg7rescue")

# set NA qval to 1
for (k in paste0(c(kd_cond,rescue_cond),"_qval")){
	tab[which(is.na(tab[,k])),k] = 1
}

# set "negative" genes qval to 1
for (k in kd_cond){
	b=paste0(k,"_b")
	qval=paste0(k,"_qval")
	tab[which(tab[,b]<minFC),qval] = 1
}
for (k in rescue_cond){
	b=paste0(k,"_b")
	qval=paste0(k,"_qval")
	tab[which(tab[,b]>-minFC),qval] = 1
}

# meta analysis
tab$UPF1 = apply(tab[,c("upf1_kd_qval","upf1_rescue_qval")], 1, function(x) psumunif(x,2))
tab$SMG6 = apply(tab[,c("smg6_kd_qval","smg6_rescue_qval")], 1, function(x) psumunif(x,2))
tab$SMG7 = apply(tab[,c("smg7_kd_qval","smg7_rescue_qval")], 1, function(x) psumunif(x,2))
tab$dKD_SMG6 = apply(tab[,c("dKD_qval","dKD_smg6rescue_qval")], 1, function(x) psumunif(x,2))
tab$dKD_SMG7 = apply(tab[,c("dKD_qval","dKD_smg7rescue_qval")], 1, function(x) psumunif(x,2))

tab$meta_SMG6=apply(tab[,c("SMG6","dKD_SMG6")], 1, function(x) Fisher.test(x))
tab$meta_SMG7=apply(tab[,c("SMG7","dKD_SMG7")], 1, function(x) Fisher.test(x))
tab$meta_SMGs=apply(tab[,c("meta_SMG6","meta_SMG7")], 1, function(x) Fisher.test(x))

tab$meta_meta=apply(tab[,c("meta_SMGs","UPF1")], 1, function(x) psumunif(x,2))

tab=tab[order(tab$meta_meta,decreasing=F),]
write.table(file="meta",tab,quote=F,row.names=F)
