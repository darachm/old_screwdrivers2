
```{r,echo=F}
#by Darach Miller
```

```{r, cache=T}
rawcounts <- read.delim("exp127counts150821_19:25:50.txt",row.names=1)
rownames(rawcounts)[!grepl("Y.{5}[WC].*",rownames(rawcounts))]
rawcounts <- rawcounts[grepl("Y.{5}[WC].*",rownames(rawcounts)),]

index <- read.csv("exp127index.csv",as.is=T)
index$time <- index$time - min(index$time)
rownames(index) <- index$sample
subrawcounts <- rawcounts[,as.logical(apply(
	sapply(paste0("Sample",index$sample,"_[(UP)|(DOWN)]"),
		grepl,names(rawcounts)),1,sum))]
sc <- subrawcounts

dubious <- read.csv("~/lab/data/otherdatasets/sgd_dubious_orfz_150831.csv",header=F)
dubiousorfz <- c("YDL227C",as.character(dubious$V2)) 
```

```{r,cache=T,warning=F}
ggplot(melt(sc))+aes(x=log(value),col=variable)+theme_bw()+
	geom_density()+guides(color=F)
ggplot(data.frame(sample=names(sc),total_counts=apply(sc,2,sum)))+
	theme_bw()+aes(x=sample,y=total_counts,label=sample)+
	geom_text(size=3)+theme(axis.text.x=element_blank())
ggplot(data.frame(sample=names(sc),median_count=apply(sc,2,median)))+
	theme_bw()+aes(x=sample,y=median_count,label=sample)+
	geom_text(size=3)+theme(axis.text.x=element_blank())
```

```{r,cache=T}
require(reshape2)
msc <- melt(data.frame(syst=rownames(sc),sc),id.vars="syst")
msc$sample <- sub("Sample(\\d+)_.+","\\1",msc$variable)
msc$tag <- factor(sub("Sample\\d+_([(UP)(DOWN)])","\\1",msc$variable))
msc$replicate <- factor(index[msc$sample,"replicate"])
msc$pulse <- factor(index[msc$sample,"pulse"])
msc$time <- as.numeric(as.character(index[msc$sample,"time"]))
dim(msc)
```

Let's use this to look at some QC before normalization.

```{r,cache=T}
dmsc <- dcast(msc,sample+replicate+pulse+time+syst~tag)
rawcorzp <- c()
rawcorzs <- c()
for (samz in unique(dmsc$sample)) {
	sdmsc <- subset(dmsc,sample==samz)
	rawcorzp[samz] <- signif(cor(sdmsc$UP,sdmsc$DOWN,
					method="pearson"),digits=2)
	rawcorzs[samz] <- signif(cor(sdmsc$UP,sdmsc$DOWN,
					method="spearman"),digits=2)
}

ggplot(data.frame(sample=names(rawcorzp),pearson=rawcorzp,spearman=rawcorzs))+
	aes(x=pearson,y=spearman,label=sample)+
	theme_bw()+geom_text(angle=-30)
subset(index,sample%in%names(which(rawcorzp<.4)))
subset(index,sample%in%names(which(rawcorzs<.4)))
```

Something looks up with sample # 5.

```{r,cache=T,fig.height=30}
dmsc$corzp <- paste("Pearson",rawcorzp[dmsc$sample])
dmsc$corzs <- paste("Spearman",rawcorzs[dmsc$sample])
ggplot(dmsc)+theme_bw()+aes(x=log(DOWN),y=log(UP))+
	geom_point(cex=.05,alpha=0.5)+
	geom_text(data=subset(dmsc,syst=="YKR039W"),x=4,y=9,aes(label=corzp))+
	geom_text(data=subset(dmsc,syst=="YKR039W"),x=4,y=8,aes(label=corzs))+
	facet_wrap(sample~replicate+pulse,ncol=4)
```

Which ones do we drop?

```{r,cache=F}
msc <- subset(msc,!(sample%in%c(73,80,85,97,109)))
msc <- msc[!(msc$variable%in%paste0("Sample",
	subset(index,(pulse==F&replicate=="b"&time>20))$sample,"_",c("UP","DOWN"))),]
```

I drop all of the sample because without the UP tag, the counts are going to be
much lower for that timepoint because of a lack of sum. Can't have that.
`msc` is a `data.frame` with a row for each observation in this experiment.
Turn it back into a dataframe and norm it.

Let's check out what dubious are doing.

```{r,cache=F}
dubz <- dcast(subset(msc,syst%in%dubiousorfz),syst~variable)
plot(apply(dubz[,-1],2,function(x){mean(na.omit(x))}),apply(dubz[,-1],2,function(x){median(na.omit(x))}),xlab="Dubious mean",ylab="Dubious median")
```

And fuck it, let's just use the median.

```{r,cache=F}
dmsc <- dcast(msc,syst~sample+tag)
rownames(dmsc) <- dmsc$syst
dmsc <- dmsc[,-1]

subdmsc <- dcast(subset(msc,time<(0.00+3*log(2)/0.12)),syst~sample+tag)
rownames(subdmsc) <- subdmsc$syst
subdmsc <- subdmsc[,-1]

#require(edgeR)
#normdmsc <- calcNormFactors(DGEList(counts=dmsc))$samples$norm.factor 
dubz <- dcast(subset(msc,syst%in%dubiousorfz),syst~variable)
dmscmedz <- apply(dubz[,-1],2,function(x){median(na.omit(x))})
normdmsc <- 1/(dmscmedz/mean(dmscmedz))
hist(normdmsc,100)
ndmsc <- dmsc
for (i in 1:ncol(ndmsc)) {
    ndmsc[,i] <- ndmsc[,i] * normdmsc[i]
}

#normsubdmsc <- calcNormFactors(DGEList(counts=subdmsc))$samples$norm.factor 
subdubz <- dcast(subset(msc,(time<(0.00+3*log(2)/0.12)&syst%in%dubiousorfz),syst~variable))
subdmscmedz <- apply(subdubz[,-1],2,function(x){median(na.omit(x))})
normsubdmsc <- 1/(subdmscmedz/mean(subdmscmedz))
hist(normsubdmsc,100)
nsubdmsc <- subdmsc
for (i in 1:ncol(nsubdmsc)) {
    nsubdmsc[,i] <- nsubdmsc[,i] * normsubdmsc[i]
}

ggplot(data.frame(sample=names(dmsc),pre_norm=apply(dmsc,2,sum),
		normed=apply(ndmsc,2,sum)))+
	theme_bw()+aes(x=pre_norm,y=normed,label=sample)+geom_text(size=2,angle=30)+
	ggtitle("comparing means")
ggplot(data.frame(sample=names(dmsc),pre_norm=apply(dmsc,2,median),
		normed=apply(ndmsc,2,median)))+
	theme_bw()+aes(x=pre_norm,y=normed,label=sample)+geom_text(size=2,angle=30)+
	ggtitle("comparing medians")
ggplot(data.frame(sample=names(subdmsc),pre_norm=apply(subdmsc,2,sum),
		normed=apply(nsubdmsc,2,sum)))+
	theme_bw()+aes(x=pre_norm,y=normed,label=sample)+geom_text(size=2,angle=30)+
	ggtitle("comparing means")
ggplot(data.frame(sample=names(subdmsc),pre_norm=apply(subdmsc,2,median),
		normed=apply(nsubdmsc,2,median)))+
	theme_bw()+aes(x=pre_norm,y=normed,label=sample)+geom_text(size=2,angle=30)+
	ggtitle("comparing medians")

```


```{r,cache=F}
mnc <- melt(data.frame(syst=rownames(ndmsc),ndmsc,check.names=F),id.vars="syst")
mnc[,c("sample","tag")] <- colsplit(mnc[,"variable"],"_",
	c("sample","tag"))
mnc$sample <- factor(mnc$sample)
mnc$tag <- factor(mnc$tag)
mnc$replicate <- factor(index[as.character(mnc$sample),"replicate"])
mnc$pulse <- factor(index[as.character(mnc$sample),"pulse"])
mnc$time <- as.numeric(as.character(index[as.character(mnc$sample),"time"]))
head(mnc)

ggplot(mnc)+theme_bw()+aes(x=log10(value+1))+geom_histogram(binwidth=0.05)+
	geom_vline(xintercept=log10(10))

ggplot(subset(mnc,syst=="YMR199W"&value>=0))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,0650))
ggplot(subset(mnc,syst=="YMR199W"&value>50))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,0650))
ggplot(subset(mnc,syst=="YKR039W"&value>=0))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,1250))
ggplot(subset(mnc,syst=="YKR039W"&value>50))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,1250))
ggplot(subset(mnc,syst=="YOR348C"&value>=0))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,0750))
ggplot(subset(mnc,syst=="YOR348C"&value>50))+theme_bw()+
	aes(x=time,y=value,col=pulse,group=replicate:pulse:tag)+
	geom_line()+ylim(c(0,0750))

ggplot(subset(mnc,value>010))+theme_bw()+aes(x=log(value))+
	geom_histogram(binwidth=0.1)
```

To illustrate again:

```{r,cache=F}
dmnc <- dcast(subset(mnc,value>=0),syst~sample+tag)
rownames(dmnc) <- dmnc$syst
dmnc <- dmnc[,-1]
require(limma)
dev.new();pmdmnc <- plotMDS(dmnc);dev.off()
```
```{r,cache=F}
tmp <- data.frame(variable=names(pmdmnc$x),ex=pmdmnc$x,why=pmdmnc$y,
	colsplit(names(pmdmnc$x),"_",c("sample","tag")))
tmp$time <- index[as.character(tmp$sample),"time"]
tmp$rep <- index[as.character(tmp$sample),"replicate"]
tmp$pulse <- index[as.character(tmp$sample),"pulse"]
ggplot(tmp)+theme_bw()+aes(x=ex,y=why,label=time)+
	geom_line(aes(group=factor(sample),linetype=factor(rep):factor(pulse)))+
	geom_point(aes(color=tag))+
	geom_text()+
	ggtitle("Taking all")
```

```{r,cache=F}
dmnc <- dcast(subset(mnc,value>10),syst~sample+tag)
rownames(dmnc) <- dmnc$syst
dmnc <- dmnc[,-1]
dev.new()
pmdmnc <- plotMDS(dmnc)
dev.off()
```
```{r,cache=F}
tmp <- data.frame(variable=names(pmdmnc$x),ex=pmdmnc$x,why=pmdmnc$y,
	colsplit(names(pmdmnc$x),"_",c("sample","tag")))
tmp$time <- index[as.character(tmp$sample),"time"]
tmp$rep <- index[as.character(tmp$sample),"replicate"]
tmp$pulse <- index[as.character(tmp$sample),"pulse"]
ggplot(tmp)+theme_bw()+aes(x=ex,y=why,label=time)+
	geom_line(aes(group=factor(sample),linetype=factor(rep):factor(pulse)))+
	geom_point(aes(color=tag))+
	geom_text()+
	ggtitle("Taking all more than 10")
```

```{r,cache=F}
dmnc <- dcast(subset(mnc,value>50),syst~sample+tag)
rownames(dmnc) <- dmnc$syst
dmnc <- dmnc[,-1]
dev.new()
pmdmnc <- plotMDS(dmnc)
dev.off()
```
```{r,cache=F}
tmp <- data.frame(variable=names(pmdmnc$x),ex=pmdmnc$x,why=pmdmnc$y,
	colsplit(names(pmdmnc$x),"_",c("sample","tag")))
tmp$time <- index[as.character(tmp$sample),"time"]
tmp$rep <- index[as.character(tmp$sample),"replicate"]
tmp$pulse <- index[as.character(tmp$sample),"pulse"]
ggplot(tmp)+theme_bw()+aes(x=ex,y=why,label=time)+
	geom_line(aes(group=factor(sample),linetype=factor(rep):factor(pulse)))+
	geom_point(aes(color=tag))+
	geom_text()+
	ggtitle("Taking all more than 50")
```

It looks like dropping everything below 10 is okay. There's still a difference
between up and down tags, but the biological variation is in both the axis.

```{r,cache=F}
tmpmnc <- subset(mnc,value>10)
dmnc <- dcast(tmpmnc,syst~sample+tag)
rownames(dmnc) <- dmnc$syst
dmnc <- dmnc[,-1]

write.table(tmpmnc,file="qcd_counts_melty_150902.tab",quote=F,sep="\t")



submnc <- melt(data.frame(syst=rownames(nsubdmsc),nsubdmsc,check.names=F),
	id.vars="syst")
submnc[,c("sample","tag")] <- colsplit(submnc[,"variable"],"_",
	c("sample","tag"))
submnc$sample <- factor(submnc$sample)
submnc$tag <- factor(submnc$tag)
submnc$replicate <- factor(index[as.character(submnc$sample),"replicate"])
submnc$pulse <- factor(index[as.character(submnc$sample),"pulse"])
submnc$time <- as.numeric(as.character(index[as.character(submnc$sample),"time"]))
tmpsubmnc <- subset(submnc,value>10)
write.table(tmpsubmnc,file="qcd_counts_firsteight_150902.tab",quote=F,sep="\t")


```


