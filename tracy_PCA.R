library(ggplot2)
library(vegan)
library(pheatmap)

# change the line below to reflect location of PopGenClass directory
setwd("~/Dropbox/popgen2020/PopGenClass/")
mice=read.csv("Tracy_mouse_data.csv")
site=mice$Site

# reading IBS matrix
ibs=as.matrix(read.table("angsd-2-5-20.ibsMat"))
dimnames(ibs)=list(mice$ID,mice$ID)
# heatmap of the matrix, with hierarchical clustering
pheatmap(1-ibs)

# unconstrained ordination (PCoA)
pp=capscale(ibs~1)

# eigenvalues:
plot(pp$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp$CA$eig/sum(pp$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues

axes2plot=c(1,2)
scores=data.frame(pp$CA$u[,axes2plot])

# plotting in ggplot:
ggplot(scores,aes(scores[,1],scores[,2],color=site)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])+coord_equal()

# testing significance of by-site divergence
adonis(ibs~site)

# plotting vegan way
library(RColorBrewer)
myColors = brewer.pal(length(levels(site)),"Dark2")
names(myColors) = levels(site)
plot(scores[,1:2],pch=16,col=myColors[site],asp=1)
ordispider(scores[,1:2],group=site,col=myColors)
ordiellipse(scores[,1:2],group=site,col=myColors,draw="polygon",label=T)

# plotting constrained ordination (by site)
pp1=capscale(ibs~site)
axes2plot=c(1,2)
scores=data.frame(pp1$CCA$wa[,axes2plot])
ggplot() + geom_point(data=scores,aes(scores[,1],scores[,2],color=site),alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])

# ordination after taking AWAY effect of site (imagine site was a nuisance covariate)
scores=data.frame(pp1$CA$u[,axes2plot]) # after taking out effect of site
ggplot(scores,aes(scores[,1],scores[,2],color=site)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])


