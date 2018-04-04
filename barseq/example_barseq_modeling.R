
```{r,echo=F}
#by Darach Miller
```

Reading in the qc'd and normalized observed counts of each mutant at each time
in each condition. I don't care how these got here, but let's give it a go.

This table of data also already has the indexing bits added, so things like
time or condition.

```{r,cache=T,cache.lazy=F}
tmp <- read.csv("systindex150825.csv",as.is=T)
systindex <- tmp[,2]
names(systindex) <- tmp[,1]
systindex[systindex==""] <- names(systindex[systindex==""])

obz <- read.delim(file="qcd_counts_melty_150830.tab")
head(obz)
hist(log10(obz$value),20)
obz$replicate <- factor(obz$replicate)
obz$pulse <- factor(obz$pulse)
obz$tag <- factor(obz$tag)
obz$common <- systindex[as.character(obz$syst)]

subobz <- read.delim(file="qcd_counts_firsteight_150902.tab")
head(subobz)
hist(log10(subobz$value),20)
subobz$replicate <- factor(subobz$replicate)
subobz$pulse <- factor(subobz$pulse)
subobz$tag <- factor(subobz$tag)
subobz$common <- systindex[as.character(subobz$syst)]

```

Remember that problem with chemostat #5? Last timepoint of the night, set the
dilution rate to 1 instead of ~10. Yep, let's drop that vessel after that
error.

```{r,cache=T,cache.lazy=F}
obz <- subset(obz,!(pulse==F&replicate=="b"&time>20))
subobz <- subset(subobz,!(pulse==F&replicate=="b"&time>20))
```

```{r,cache=T,cache.lazy=F}
datar <- dcast(obz,formula=pulse+replicate+syst+time~tag,value.var="value")
datar$sumz <- apply(datar,1,function(x){sum(as.numeric(x[5:6]),na.rm=T)})

mutantdatar <- dcast(subset(datar),syst~pulse+replicate+time,value.var="sumz")

mdatar <- melt(datar[,-(5:6)],value.name="sumz",
	id.vars=c("pulse","replicate","syst","time"))
mdatar$common <- systindex[as.character(mdatar$syst)]
names(mdatar)[6] <- "sumz"
mdatar$replicate <- factor(mdatar$replicate)
mdatar$pulse <- factor(mdatar$pulse)

subdatar <- dcast(subobz,formula=pulse+replicate+syst+time~tag,value.var="value")
subdatar$sumz <- apply(subdatar,1,function(x){sum(as.numeric(x[5:6]),na.rm=T)})

mutantsubdatar <- dcast(subset(subdatar),syst~pulse+replicate+time,value.var="sumz")

msubdatar <- melt(subdatar[,-(5:6)],value.name="sumz",
	id.vars=c("pulse","replicate","syst","time"))
msubdatar$common <- systindex[as.character(msubdatar$syst)]
names(msubdatar)[6] <- "sumz"
msubdatar$replicate <- factor(msubdatar$replicate)
msubdatar$pulse <- factor(msubdatar$pulse)

```

```{r,cache=T,cache.lazy=F,fig.height=10}
ggplot(mdatar)+theme_bw()+aes(x=log10(sumz))+
	geom_histogram(binwidth=0.1)+facet_grid(time~pulse+replicate)
```

```{r,cache=T,cache.lazy=F,fig.height=10}
ggplot(mdatar)+theme_bw()+aes(x=(sumz))+
	geom_histogram(binwidth=50)+facet_grid(time~pulse+replicate)+
	xlim(0,4e3)
```

Okay, now for modeling. I added glutamine to the `pulse==T` samples at 1330 to 
1333, so essentially t+1.104 hours.

Here's a simple log linear model, with a term for a linear effect after adding
glutamine or not.

```{r,cache=T,cache.lazy=F,warning=F,error=F}
lmz <- apply(mutantdatar,1,function(x){
	syst <- unlist(x["syst"])
	x <- unlist(x[-1])
	xz <- data.frame(
			colsplit(names(x),"_",c("pulse","rep","time")),
			logcountz=log(as.numeric(unlist(x))))
	tmp <- try(lm(logcountz~time+time:(pulse&(time>1.1))+pulse:rep,data=xz))
	if (is(tmp)[1]=="lm") {
		if (length(summary(tmp)$coefficients)<12){return(NA)}
		return(c(as.character(syst),
			unlist(dcast(melt(summary(tmp)$coefficient),.~X1+X2)[-1])))
	} else {
		return(NA)
	}
})
lmz <- data.frame(t(data.frame(sapply(lmz,function(x){
	return(c(
		x[1],
		x["(Intercept)_Estimate"],
		x["(Intercept)_Pr(>|t|)"],
		x["(Intercept)_t value"],
		x["(Intercept)_Std. Error"],
		x["pulseFALSE:repa_Estimate"],
		x["pulseFALSE:repa_Pr(>|t|)"],
		x["pulseFALSE:repa_t value"],
		x["pulseFALSE:repa_Std. Error"],
		x["pulseFALSE:repb_Estimate"],
		x["pulseFALSE:repb_Pr(>|t|)"],
		x["pulseFALSE:repb_t value"],
		x["pulseFALSE:repb_Std. Error"],
		x["pulseTRUE:repa_Estimate"],
		x["pulseTRUE:repa_Pr(>|t|)"],
		x["pulseTRUE:repa_t value"],
		x["pulseTRUE:repa_Std. Error"],
		x["time_Estimate"],
		x["time_Pr(>|t|)"],
		x["time_t value"],
		x["time_Std. Error"],
		x["time:pulse & (time > 1.1)TRUE_Estimate"],
		x["time:pulse & (time > 1.1)TRUE_Pr(>|t|)"],
		x["time:pulse & (time > 1.1)TRUE_t value"],
		x["time:pulse & (time > 1.1)TRUE_Std. Error"]
	))
}))))
lmz <- lmz[!is.na(lmz[,1]),]
lmz <- data.frame(apply(lmz,2,as.character))
rownames(lmz)<-as.character(lmz[,1])
lmz <- data.frame(apply(lmz[,-1],2,as.numeric),row.names=rownames(lmz))
colnames(lmz) <- c(
	"intercept_estimate",
	"intercept_pval",
	"intercept_tval",
	"intercept_stder",
	"pulse.repa_estimate",
	"pulse.repa_pval",
	"pulse.repa_tval",
	"pulse.repa_stder",
	"pulse.repb_estimate",
	"pulse.repb_pval",
	"pulse.repb_tval",
	"pulse.repb_stder",
	"pulseT.repa_estimate",
	"pulseT.repa_pval",
	"pulseT.repa_tval",
	"pulseT.repa_stder",
	"time_estimate",
	"time_pval",
	"time_tval",
	"time_stder",
	"time.pulse_estimate",
	"time.pulse_pval",
	"time.pulse_tval",
	"time.pulse_stder"
)
apply(lmz,2,function(x){range(na.omit(x))})

sublmz <- apply(mutantsubdatar,1,function(x){
	syst <- unlist(x["syst"])
	x <- unlist(x[-1])
	xz <- data.frame(
			colsplit(names(x),"_",c("pulse","rep","time")),
			logcountz=log(as.numeric(unlist(x))))
	tmp <- try(lm(logcountz~time+time:(pulse&(time>1.1))+pulse:rep,data=xz))
	if (is(tmp)[1]=="lm") {
		if (length(summary(tmp)$coefficients)<12){return(NA)}
		return(c(as.character(syst),
			unlist(dcast(melt(summary(tmp)$coefficient),.~X1+X2)[-1])))
	} else {
		return(NA)
	}
})
sublmz <- data.frame(t(data.frame(sapply(sublmz,function(x){
	return(c(
		x[1],
		x["(Intercept)_Estimate"],
		x["(Intercept)_Pr(>|t|)"],
		x["(Intercept)_t value"],
		x["(Intercept)_Std. Error"],
		x["pulseFALSE:repa_Estimate"],
		x["pulseFALSE:repa_Pr(>|t|)"],
		x["pulseFALSE:repa_t value"],
		x["pulseFALSE:repa_Std. Error"],
		x["pulseFALSE:repb_Estimate"],
		x["pulseFALSE:repb_Pr(>|t|)"],
		x["pulseFALSE:repb_t value"],
		x["pulseFALSE:repb_Std. Error"],
		x["pulseTRUE:repa_Estimate"],
		x["pulseTRUE:repa_Pr(>|t|)"],
		x["pulseTRUE:repa_t value"],
		x["pulseTRUE:repa_Std. Error"],
		x["time_Estimate"],
		x["time_Pr(>|t|)"],
		x["time_t value"],
		x["time_Std. Error"],
		x["time:pulse & (time > 1.1)TRUE_Estimate"],
		x["time:pulse & (time > 1.1)TRUE_Pr(>|t|)"],
		x["time:pulse & (time > 1.1)TRUE_t value"],
		x["time:pulse & (time > 1.1)TRUE_Std. Error"]
	))
}))))
sublmz <- sublmz[!is.na(sublmz[,1]),]
sublmz <- data.frame(apply(sublmz,2,as.character))
rownames(sublmz)<-as.character(sublmz[,1])
sublmz <- data.frame(apply(sublmz[,-1],2,as.numeric),row.names=rownames(sublmz))
colnames(sublmz) <- c(
	"intercept_estimate",
	"intercept_pval",
	"intercept_tval",
	"intercept_stder",
	"pulse.repa_estimate",
	"pulse.repa_pval",
	"pulse.repa_tval",
	"pulse.repa_stder",
	"pulse.repb_estimate",
	"pulse.repb_pval",
	"pulse.repb_tval",
	"pulse.repb_stder",
	"pulseT.repa_estimate",
	"pulseT.repa_pval",
	"pulseT.repa_tval",
	"pulseT.repa_stder",
	"time_estimate",
	"time_pval",
	"time_tval",
	"time_stder",
	"time.pulse_estimate",
	"time.pulse_pval",
	"time.pulse_tval",
	"time.pulse_stder"
)
apply(sublmz,2,function(x){range(na.omit(x))})

```

```{r,cache=T,cache.lazy=F}
ggplot(lmz)+aes(x=time_tval)+geom_histogram(binwidth=0.4)
```

Okay, so that looks really funny. 

What about how a chemostat doesn't select for maximal growth rate, it selects
for staying in the vessel with a constant dilution rate. Population size is not
constant unless you believe yeasts are elbowing each other out to the dilution
tube, so maybe increased abundance of some mutants can elbow others out of the
sequencing normalization. Spike-ins would be nice right now.

What about using dubious ORFs as a spike in? I queried SGD for all dubious ORFs
and added on the HO locus. Let's see what those do.

```{r,cache=T,cache.lazy=F}
dubious <- read.csv("~/lab/data/otherdatasets/sgd_dubious_orfz_150831.csv",header=F)
dubiousorfz <- c("YDL227C",as.character(dubious$V2))
ggplot(subset(lmz,rownames(lmz)%in%dubiousorfz))+
	geom_histogram(binwidth=0.005)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs, all")+aes(x=time_estimate)
ggplot(subset(lmz,rownames(lmz)%in%dubiousorfz&time_pval<0.05))+
	geom_histogram(binwidth=0.005)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs, pval < 0.05")+aes(x=time_estimate)
basaldilution <- median(subset(lmz,rownames(lmz)%in%dubiousorfz&time_pval<0.05)$time_estimate)
subbasaldilution <- median(subset(sublmz,rownames(sublmz)%in%dubiousorfz&time_pval<0.05)$time_estimate)
summary(subset(lmz,rownames(lmz)%in%dubiousorfz&time_pval<0.05)$time_estimate)
```

So we use dubious orfs to find out how much the faster growing stuff is
probably diluting everything else.
There's also the issue of the pulse.

```{r,cache=T,cache.lazy=F}
ggplot(subset(lmz,rownames(lmz)%in%dubiousorfz))+
	geom_histogram(binwidth=0.0010)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs, all")+aes(x=time.pulse_estimate)
ggplot(subset(lmz,rownames(lmz)%in%dubiousorfz&time.pulse_pval<0.05))+
	geom_histogram(binwidth=0.0010)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs, pval < 0.05")+aes(x=time.pulse_estimate)
basalpulse <- median(subset(lmz,
	rownames(lmz)%in%dubiousorfz&time.pulse_pval<0.05)$time.pulse_estimate)
subbasalpulse <- median(subset(sublmz,
	rownames(sublmz)%in%dubiousorfz&time.pulse_pval<0.05)$time.pulse_estimate)
summary(subset(lmz,rownames(lmz)%in%dubiousorfz&time.pulse_pval<0.05)$time.pulse_estimate)
```

Most things grows faster with a little bit more glutamine.

Below, we re-do the log-linear modeling but we add back the median of the
estimated rates for the dubious orfs that fit the model at pval < 0.05, and
we add back some for the pulse.

My null hypothesis is the dubious orfs.

```{r,cache=T,cache.lazy=F,warning=F,error=F}
adjlmz <- apply(mutantdatar,1,function(x){
	syst <- unlist(x["syst"])
	x <- unlist(x[-1])
	xz <- data.frame(
			colsplit(names(x),"_",c("pulse","rep","time")),
			logcountz=log(as.numeric(unlist(x))))
	#### FANCY BIT ####
	xz$logcountz <- xz$logcountz + -basaldilution*xz$time
	xz$logcountz[xz$pulse&xz$time>1.1] <- xz$logcountz[xz$pulse&xz$time>1.1] + 
		-basalpulse*xz$time[xz$pulse&xz$time>1.1]
	#### FANCY BIT ####
	tmp <- try(lm(logcountz~time+time:(pulse&(time>1.1))+pulse:rep,data=xz))
	if (is(tmp)[1]=="lm") {
		if (length(summary(tmp)$coefficients)<12){return(NA)}
		return(c(as.character(syst),
			unlist(dcast(melt(summary(tmp)$coefficient),.~X1+X2)[-1])))
	} else {
		return(NA)
	}
})
adjlmz <- data.frame(t(data.frame(sapply(adjlmz,function(x){
	return(c(
		x[1],
		x["(Intercept)_Estimate"],
		x["(Intercept)_Pr(>|t|)"],
		x["(Intercept)_t value"],
		x["(Intercept)_Std. Error"],
		x["pulseFALSE:repa_Estimate"],
		x["pulseFALSE:repa_Pr(>|t|)"],
		x["pulseFALSE:repa_t value"],
		x["pulseFALSE:repa_Std. Error"],
		x["pulseFALSE:repb_Estimate"],
		x["pulseFALSE:repb_Pr(>|t|)"],
		x["pulseFALSE:repb_t value"],
		x["pulseFALSE:repb_Std. Error"],
		x["pulseTRUE:repa_Estimate"],
		x["pulseTRUE:repa_Pr(>|t|)"],
		x["pulseTRUE:repa_t value"],
		x["pulseTRUE:repa_Std. Error"],
		x["time_Estimate"],
		x["time_Pr(>|t|)"],
		x["time_t value"],
		x["time_Std. Error"],
		x["time:pulse & (time > 1.1)TRUE_Estimate"],
		x["time:pulse & (time > 1.1)TRUE_Pr(>|t|)"],
		x["time:pulse & (time > 1.1)TRUE_t value"],
		x["time:pulse & (time > 1.1)TRUE_Std. Error"]
	))
}))))
adjlmz <- adjlmz[!is.na(adjlmz[,1]),]
adjlmz <- data.frame(apply(adjlmz,2,as.character))
rownames(adjlmz)<-as.character(adjlmz[,1])
adjlmz <- data.frame(apply(adjlmz[,-1],2,as.numeric),row.names=rownames(adjlmz))
colnames(adjlmz) <- c(
	"intercept_estimate",
	"intercept_pval",
	"intercept_tval",
	"intercept_stder",
	"pulse.repa_estimate",
	"pulse.repa_pval",
	"pulse.repa_tval",
	"pulse.repa_stder",
	"pulse.repb_estimate",
	"pulse.repb_pval",
	"pulse.repb_tval",
	"pulse.repb_stder",
	"pulseT.repa_estimate",
	"pulseT.repa_pval",
	"pulseT.repa_tval",
	"pulseT.repa_stder",
	"time_estimate",
	"time_pval",
	"time_tval",
	"time_stder",
	"time.pulse_estimate",
	"time.pulse_pval",
	"time.pulse_tval",
	"time.pulse_stder"
)
apply(adjlmz,2,function(x){range(na.omit(x))})

subadjlmz <- apply(mutantsubdatar,1,function(x){
	syst <- unlist(x["syst"])
	x <- unlist(x[-1])
	xz <- data.frame(
			colsplit(names(x),"_",c("pulse","rep","time")),
			logcountz=log(as.numeric(unlist(x))))
	#### FANCY BIT ####
	xz$logcountz <- xz$logcountz + -basaldilution*xz$time
	xz$logcountz[xz$pulse&xz$time>1.1] <- xz$logcountz[xz$pulse&xz$time>1.1] + 
		-basalpulse*xz$time[xz$pulse&xz$time>1.1]
	#### FANCY BIT ####
	tmp <- try(lm(logcountz~time+time:(pulse&(time>1.1))+pulse:rep,data=xz))
	if (is(tmp)[1]=="lm") {
		if (length(summary(tmp)$coefficients)<12){return(NA)}
		return(c(as.character(syst),
			unlist(dcast(melt(summary(tmp)$coefficient),.~X1+X2)[-1])))
	} else {
		return(NA)
	}
})
subadjlmz <- data.frame(t(data.frame(sapply(subadjlmz,function(x){
	return(c(
		x[1],
		x["(Intercept)_Estimate"],
		x["(Intercept)_Pr(>|t|)"],
		x["(Intercept)_t value"],
		x["(Intercept)_Std. Error"],
		x["pulseFALSE:repa_Estimate"],
		x["pulseFALSE:repa_Pr(>|t|)"],
		x["pulseFALSE:repa_t value"],
		x["pulseFALSE:repa_Std. Error"],
		x["pulseFALSE:repb_Estimate"],
		x["pulseFALSE:repb_Pr(>|t|)"],
		x["pulseFALSE:repb_t value"],
		x["pulseFALSE:repb_Std. Error"],
		x["pulseTRUE:repa_Estimate"],
		x["pulseTRUE:repa_Pr(>|t|)"],
		x["pulseTRUE:repa_t value"],
		x["pulseTRUE:repa_Std. Error"],
		x["time_Estimate"],
		x["time_Pr(>|t|)"],
		x["time_t value"],
		x["time_Std. Error"],
		x["time:pulse & (time > 1.1)TRUE_Estimate"],
		x["time:pulse & (time > 1.1)TRUE_Pr(>|t|)"],
		x["time:pulse & (time > 1.1)TRUE_t value"],
		x["time:pulse & (time > 1.1)TRUE_Std. Error"]
	))
}))))
subadjlmz <- subadjlmz[!is.na(subadjlmz[,1]),]
subadjlmz <- data.frame(apply(subadjlmz,2,as.character))
rownames(subadjlmz)<-as.character(subadjlmz[,1])
subadjlmz <- data.frame(apply(subadjlmz[,-1],2,as.numeric),row.names=rownames(subadjlmz))
colnames(subadjlmz) <- c(
	"intercept_estimate",
	"intercept_pval",
	"intercept_tval",
	"intercept_stder",
	"pulse.repa_estimate",
	"pulse.repa_pval",
	"pulse.repa_tval",
	"pulse.repa_stder",
	"pulse.repb_estimate",
	"pulse.repb_pval",
	"pulse.repb_tval",
	"pulse.repb_stder",
	"pulseT.repa_estimate",
	"pulseT.repa_pval",
	"pulseT.repa_tval",
	"pulseT.repa_stder",
	"time_estimate",
	"time_pval",
	"time_tval",
	"time_stder",
	"time.pulse_estimate",
	"time.pulse_pval",
	"time.pulse_tval",
	"time.pulse_stder"
)
apply(subadjlmz,2,function(x){range(na.omit(x))})

```

```{r,cache=T,cache.lazy=F}
ggplot(adjlmz)+aes(x=time_tval)+geom_histogram(binwidth=0.4)
ggplot(adjlmz)+aes(x=time.pulse_tval)+geom_histogram(binwidth=0.1)
ggplot(subset(adjlmz,rownames(adjlmz)%in%dubiousorfz))+
	geom_histogram(binwidth=0.001)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs")+aes(x=time_estimate)
ggplot(subset(adjlmz,rownames(adjlmz)%in%dubiousorfz))+
	geom_histogram(binwidth=0.001)+xlim(-0.1,0.05)+
	ggtitle("Just the dubious orfs")+aes(x=time.pulse_estimate)

ggplot(subset(lmz,time_pval<0.05))+aes(x=time_estimate)+
	geom_histogram(binwidth=0.001)+
	xlim(-0.3,0.10)+ylim(0,210)+
	ggtitle("Pre-adjust, All fitness coefficients, pval < 0.05")
ggplot(subset(adjlmz,time_pval<0.05))+aes(x=time_estimate)+
	geom_histogram(binwidth=0.001)+
	xlim(-0.3,0.10)+ylim(0,210)+
	ggtitle("Adjusted, All fitness coefficients, pval < 0.05")
ggplot(subset(adjlmz,time.pulse_pval<0.05))+aes(x=time.pulse_estimate)+
	geom_histogram(binwidth=0.001)+
	xlim(-0.3,0.10)+ylim(0,210)+
	ggtitle("Adjusted, All pulse fitness coefficients, pval < 0.05")


```

```{r,cache=T,cache.lazy=F}
ggplot(adjlmz)+aes(x=log10(time.pulse_pval))+geom_histogram(binwidth=0.1)
ggplot(adjlmz)+aes(x=log10(time_pval))+geom_histogram(binwidth=0.1)

sigtimeadjlmz <- adjlmz[!is.na(adjlmz$time_pval),]
sigtimeadjlmz <- sigtimeadjlmz[qvalue(sigtimeadjlmz$time_pval,fdr.level=0.05)$significant,]
ggplot(sigtimeadjlmz)+aes(x=time_estimate)+geom_histogram(binwidth=0.001)
```

```{r,cache=T,cache.lazy=F,fig.height=400}
ggplot(subset(mdatar,syst%in%rownames(subset(sigtimeadjlmz,time_estimate>0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=3)+
	ggtitle("Positive time_estimate")
```
```{r,cache=T,cache.lazy=F,fig.height=200}
ggplot(subset(mdatar,syst%in%rownames(subset(sigtimeadjlmz,time_estimate<0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=3)+
	ggtitle("Negative time_estimate")
```

```{r,cache=T,cache.lazy=F}
#
sigtime.pulseadjlmz <- adjlmz[!is.na(adjlmz$time.pulse_pval),]
sigtime.pulseadjlmz <- sigtime.pulseadjlmz[qvalue(sigtime.pulseadjlmz$time.pulse_pval,fdr.level=0.05)$significant,]
ggplot(sigtime.pulseadjlmz)+aes(x=time.pulse_estimate)+geom_histogram(binwidth=0.01)
```

```{r,cache=T,cache.lazy=F,fig.height=05}
ggplot(subset(mdatar,syst%in%rownames(subset(sigtime.pulseadjlmz,time.pulse_estimate>0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=2)+
	ggtitle("Positive time.pulse_estimate")
```
```{r,cache=T,cache.lazy=F,fig.height=40}
ggplot(subset(mdatar,
		syst%in%rownames(subset(sigtime.pulseadjlmz,time.pulse_estimate<0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=2)+
	ggtitle("Negative time.pulse_estimate")
```

```{r,cache=T,cache.lazy=F}
sigtime.pulsesubadjlmz <- subadjlmz[!is.na(subadjlmz$time.pulse_pval),]
sigtime.pulsesubadjlmz <- sigtime.pulsesubadjlmz[sigtime.pulsesubadjlmz$time.pulse_pval<0.05,]
#qvalue(sigtime.pulsesubadjlmz$time.pulse_pval,fdr.level=0.05)$significant,]
ggplot(sigtime.pulsesubadjlmz)+aes(x=time.pulse_estimate)+geom_histogram(binwidth=0.01)
```

```{r,cache=T,cache.lazy=F,fig.height=10}
ggplot(subset(msubdatar,syst%in%rownames(subset(sigtime.pulsesubadjlmz,time.pulse_estimate>0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=2)+
	ggtitle("Positive time.pulse_estimate")
```
```{r,cache=T,cache.lazy=F,fig.height=10}
ggplot(subset(msubdatar,
		syst%in%rownames(subset(sigtime.pulsesubadjlmz,time.pulse_estimate<0))))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~common,scales="free_y",ncol=2)+
	ggtitle("Negative time.pulse_estimate")
```

```{r,cache=T,cache.lazy=F,fig.height=10}
ggplot(subset(subdatar,syst%in%c("YNL068C","YIL131C")))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~syst,scales="free_y")
ggplot(subset(datar,syst%in%c("YNL068C","YIL131C")))+
	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
	geom_point()+geom_line()+
	facet_wrap(~syst,scales="free_y")
```



```{r,cache=T,cache.lazy=F}
#ggplot(subset(mdatar,syst%in%
#		rownames(pulserankadjlmz)[c(1:12,(nrow(pulserankadjlmz)-11):nrow(pulserankadjlmz))]))+
#	theme_bw()+aes(x=time,y=sumz,col=pulse,group=pulse:replicate)+
#	geom_point()+geom_line()+
#	facet_wrap(~common,scales="free_y",ncol=3)+
#	ggtitle("Things that should have a fitness effect, with respect to the pulse")
```

