library(ggplot2)

tab=read.table("simus")
names(tab)=c("perc","pearson")
tab$perc=tab$perc/2
#tab2=tab
#tab2$perc=tab2$perc^2
#fit=lm(perc~pearson,data=tab2)
#fitted=sqrt(predict(fit,newdata=data.frame(pearson=c(0.8424621)),interval="prediction"))
fit=loess(perc~pearson,data=tab,span=0.5)
pred=predict(fit,newdata=data.frame(pearson=c(0.8424621)),interval="prediction",se=T,level=0.95)
p=ggplot(tab,aes(x=perc,y=pearson))
svg(height=2.8,width=5.5)
	p + 
	geom_point() + 
	geom_smooth(se=F) + 
	geom_vline(xintercept=pred$fit+2*pred$se.fit) + 
	geom_vline(xintercept=pred$fit-2*pred$se.fit) + 
	geom_vline(xintercept=pred$fit) + 
	geom_hline(yintercept=0.8424621) + 
	xlab("average fold-change variation") + 
	ylab("Pearson's r") +
	theme_bw()
	#p + geom_point() + geom_smooth() + xlab("average fold-change variation") + ylab("Pearson's r")
dev.off()
print(pred$fit)
print(pred$se.fit)
