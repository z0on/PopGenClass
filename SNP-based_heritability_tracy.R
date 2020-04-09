#install.packages("MCMCglmm")
library(MCMCglmm)
library(MASS)
library(pheatmap)
library(vegan)

mice=read.csv("~/Dropbox/popgen2020/PopGenClass/Tracy_mouse_data.csv")
row.names(mice)=mice$ID
# correcting negative latitude entries :)
mice$Lat=abs(mice$Lat)
table(mice$Site)

# removing outlying individuals (prevent constructon of a good inverted IBS matrix)
# bads=which(mice$ID %in% c("SGTm082","SGTm122"))
bads=which(mice$ID %in% c("SGTm082"))
mice=mice[-bads,]
mice$Site=factor(mice$Site,levels=unique(mice$Site))

#------- making inverted genetic similarity matrix

ibs=as.matrix(read.table("~/Dropbox/popgen2020/PopGenClass/angsd-2-5-20.ibsMat"))
ibs=ibs[-bads,-bads]
dimnames(ibs)=list(mice$ID,mice$ID)
pheatmap(1-ibs,cex=0.8)
plot(capscale(ibs~1))
# making generalized inverse of IBS matrix
gc0=1-ibs
gc=as(ginv(gc0),"dgCMatrix")
dimnames(gc)=list(unique(mice$ID),unique(mice$ID))
pheatmap(gc,cex=0.7)

# check if the matrix is "nice" (positive definite)
table(eigen(gc)$values<0)

# ---- site coordinates: making inverted covariance matrix for sites

sitexy=c()
for (s in unique(mice$Site)) {
	sitexy=data.frame(rbind(sitexy,subset(mice,Site==s)[1,]))
}
xy=sitexy[,c("Lat","Lon")]
dxy=as.matrix(dist(xy,diag=T,upper = T))
# making generalized inverse of between-site distance matrix
idxy0=1-dxy
idxy=as(ginv(idxy0),"dgCMatrix")
dimnames(idxy)=list(unique(mice$Site),unique(mice$Site))
pheatmap(idxy)
# check if the matrix is "nice" (positive definite)
table(eigen(idxy)$values<0)

#---------- fitting model and analyzing heritability

# interesting traits:
names(mice)[11:20]

# heritable ones are the ones about call frequencies, with Hz or BW, of DF in them. 

trait="note_num"
mice$Resp=mice[,trait]

# "slightly informative" prior (often used in papers)- basically we are willing to assume that each variance component explains at least some variance, and we race them against each other starting from the same prior variance (1/3d of total variance) to see which explains more.

rv=var(mice$Resp,na.rm=T)
prior1=list(
	  R=list(V=rv/3,nu=1),
	  G=list(
	    G1=list(V=rv/3,nu=1),
	    G2=list(V=rv/3,nu=1)
	    )
	)
	
# uninformative priors with parameter expansion
# for model with two random effects
prior0=list(
	   R=list(V=1,nu=0.002),
	   G=list(
	     G1=list(V=1,nu=1, alpha.mu=0,alpha.V=1000),
	     G2=list(V=1,nu=1, alpha.mu=0,alpha.V=1000)
	     )
	 )
	 
# for model with one random effect:
prior00=list(
	   R=list(V=1,nu=0.002),
	   G=list(
	     G1=list(V=1,nu=1, alpha.mu=0,alpha.V=1000)
	     )
	 )
	
#fitting model with ID and site
mm=MCMCglmm(Resp~Sex,random=~ID+Site,ginverse=list(ID=gc,Site=idxy),data=mice, prior=prior0,nitt=23000,burnin=3000,thin=20)
mm1=MCMCglmm(Resp~Sex,random=~ID+Site,ginverse=list(ID=gc,Site=idxy),data=mice, prior=prior1,nitt=23000,burnin=3000,thin=20)
mm01=MCMCglmm(Resp~Sex,random=~ID+Site,ginverse=list(ID=gc),data=mice, prior=prior0,nitt=23000,burnin=3000,thin=20)
mm00=MCMCglmm(Resp~Sex,random=~ID,ginverse=list(ID=gc),data=mice, prior=prior00,nitt=23000,burnin=3000,thin=20)
mmf=MCMCglmm(Resp~Sex+Site,random=~ID,ginverse=list(ID=gc),data=mice, prior=prior00,nitt=23000,burnin=3000,thin=20)

# basic summary: note eff.samp numbers for each variance component. They might be quite small because of autocorrelation between sampled values. If so, increase nitt and thin in the call to MCMCglmm above.
summary(mm)

# fixed effects:
colnames(mm$Sol)
# fixed effect of sexM:
alpha=mm$Sol[,"SexM"]
plot(alpha)
mean(alpha)
HPDinterval(alpha)

# variance compionents:
colnames(mm$VCV)

# heritability: proportion of variation attributable to ID (first variance component):
mmm=mmf
totalVar=apply(mmm$VCV,1,sum)
h2=mmm$VCV[,"ID"]/totalVar
plot(h2)
mean(h2)
HPDinterval(h2)

# proportion of variation due to site:
sit=mm$VCV[,"Site"]/totalVar
plot(sit)
mean(sit)
HPDinterval(sit)
