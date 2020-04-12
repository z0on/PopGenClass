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
nrow(mice)

# removing outlying individuals (look at pheatmap of inverted genetic matrix, line 40)
bads=0
#bads=which(mice$ID %in% c("SGTm082","SGTf03","SGTm029","SGTm122"))
if(bads[1]!=0) { mice=mice[-bads,] }
mice$Site=factor(mice$Site,levels=unique(mice$Site))

#------- making inverted genetic similarity matrix

covmat=as.matrix(read.table("~/Dropbox/popgen2020/PopGenClass/angsd-2-5-20.covMat"))
if(bads[1]!=0) { covmat=covmat[-bads,-bads] }
dimnames(covmat)=list(mice$ID,mice$ID)
pheatmap(covmat,cex=0.8)

# ensuring the matrix is positive definite (adjustin it just a bit with nearPD function)
if (range(eigen(covmat)$values)[1]<0) {
	covmat1=as.matrix(nearPD(covmat,ensureSymmetry=T,conv.norm.type="F")$mat)
# were adjustments quite subtle?
	plot(covmat1[upper.tri(covmat1)]~ covmat[upper.tri(covmat)])
	abline(0,1,col="red")
	covmat=covmat1
}

# making generalized inverse of covmat matrix
gci=ginv(covmat)
dimnames(gci)=list(unique(mice$ID),unique(mice$ID))
pheatmap(gci,cex=0.8)
# check if some samples are really far outside everyone else. If so, list them on line 16 and redo everyting from line 7. 

# ensuring the matrix is positive definite (just in case)
gci1=as.matrix(nearPD(gci,ensureSymmetry=T,conv.norm.type="F")$mat)
# were adjustments quite subtle?
plot(gci1[upper.tri(gci1)]~ gci[upper.tri(gci)])
gc=as(gci1,"dgCMatrix")

# ------- eigenvalue decomposition (just to plot PCA fo fun)

eig=eigen(covmat)
plot(eig$vectors,asp=1)
# pcoa based on correlations
pcoa=capscale((1-cov2cor(covmat))~1)

# how different are eigens and pcoa?
plot(procrustes(eig$vectors,pcoa))


# ---- site coordinates: making inverted covariance matrix for sites

sitexy=c()
for (s in unique(mice$Site)) {
	sitexy=data.frame(rbind(sitexy,subset(mice,Site==s)[1,]))
}
xy=sitexy[,c("Lat","Lon")]
dxy=as.matrix(dist(xy,diag=T,upper = T))
dimnames(dxy)=list(sitexy$Site,sitexy$Site)
pheatmap(dxy)
# making generalized inverse of between-site distance matrix
idxy0=1-dxy
idxy=as(ginv(idxy0),"dgCMatrix")
dimnames(idxy)=list(unique(mice$Site),unique(mice$Site))
pheatmap(idxy)
# check if the matrix is "nice" (positive definite)
table(eigen(idxy)$values<0)

#--------- relatedness (skip this)

 rel=read.table("~/Dropbox/tracy_heritability/ANGSD_related.output",sep="\t",header=T)
 relm=matrix(0,nrow=length(unique(rel$a))+1,ncol=length(unique(rel$a))+1)
 dim(relm)
 for (a in unique(rel$a)) {
	 for (b in unique(rel$b)) {
		 if (b<=a) {next}
		 relm[a,b]=relm[b,a]=rel[rel$a==a & rel$b==b,"rab"]
	 }
 }
 diag(relm)=1
 relm=relm[-bads,-bads]
 pp=pheatmap(relm)
rm=as(ginv(relm),"dgCMatrix")
dimnames(rm)=list(mice$ID,mice$ID)
pp=pheatmap(rm)
# plot(pp$tree_row,cex=0.5)
table(eigen(rm)$values<0)

#---------- fitting model and analyzing heritability

# interesting traits:
names(mice)[11:19]

# heritable ones are the ones about call frequencies, with Hz or BW, of DF in them. 

trait="call_Hz_min"
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

mm01=MCMCglmm(Resp~Sex,random=~ID+Site,ginverse=list(ID=gc),data=mice, prior=prior1,nitt=23000,burnin=3000,thin=20)

mm00=MCMCglmm(Resp~Sex,random=~ID,ginverse=list(ID=gc),data=mice, prior=prior00,nitt=23000,burnin=3000,thin=20)

mmf=MCMCglmm(Resp~Sex+Site,random=~ID,ginverse=list(ID=gc),data=mice, prior=prior00,nitt=23000,burnin=3000,thin=20)

# ------- basic summary
#note eff.samp numbers for each variance component. They might be quite small because of autocorrelation between sampled values. If so, increase nitt and thin in the call to MCMCglmm above.

model=mm1
summary(model)
plot(model)

# ------- fixed effects:
colnames(model$Sol)
# fixed effect of sexM:
alpha= model$Sol[,"SexM"]
plot(alpha)
mean(alpha)
HPDinterval(alpha)

# ------- variance compionents:
colnames(model$VCV)

# ------- heritability: proportion of variation attributable to ID (first variance component):
totalVar=apply(model$VCV,1,sum)
h2= model$VCV[,"ID"]/totalVar
plot(h2)
mean(h2)
HPDinterval(h2)

# -------- proportion of variation due to site:
sit= model$VCV[,"Site"]/totalVar
plot(sit)
mean(sit)
HPDinterval(sit)
