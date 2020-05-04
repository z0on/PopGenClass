# install.packages("vegan")
library(ggplot2)
library(vegan)
library(pheatmap)

# change the line below to reflect location of PopGenClass directory
setwd("~/Dropbox/popgen2020/PopGenClass/")

# reading metadata; all we need from this now is "site" - place where mice were caught
mice=read.csv("Tracy_mouse_data.csv")
site=mice$Site

# reading IBS matrix
ibs=as.matrix(read.table("angsd-2-5-20.ibsMat"))
dimnames(ibs)=list(mice$ID,mice$ID)
# heatmap of the matrix, with hierarchical clustering (for shits and giggles)
pheatmap(ibs)

# unconstrained ordination (PCoA)
pp=capscale(ibs~1)
# simple plot:
plot(pp)

# eigenvalues: how many good eigenvalues (exceeding broken stick expectation) are there:
plot(pp$CA$eig,ylab="proportion of variance explained",xlab="eigenvector #")
plot(pp$CA$eig/sum(pp$CA$eig),ylab="proportion of variance explained",xlab="eigenvector #")
lines(bstick(length(pp$CA$eig)),col="red") # "broken stick" model: expected random distribution of eigenvalues
# looks like we have 2 good eigenvalues

# extracting "scores" table, to plot
axes2plot=c(1,2) # which PCAs to plot
scores=data.frame(pp$CA$u[,axes2plot])

# -------- plotting in ggplot, colored by site:
quartz()
ggplot(scores,aes(scores[,1],scores[,2],color=site)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])+coord_equal()

# testing significance of by-site divergence
adonis(ibs~site)
# R2 = proportion of variance attributable to the factor. Here, "site" explains 7.58% of variance
           # Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)    
# site        5    0.7717 0.154342  2.5588 0.0758  0.001 ***
# Residuals 156    9.4094 0.060317         0.9242           
# Total     161   10.1811                  1.0000           

# -------- plotting vegan way (I really like ordispider and ordiellipse functions, btu colors are a bit of trouble)

library(RColorBrewer)
# making color scheme for our sites (will not work if there are >12 "site" levels)
myColors = brewer.pal(length(levels(site)),"Set3")
names(myColors) = unique(site)
# will work with many more "site" levels, if you have WGCNA:
# myColors  = WGCNA::labels2colors(levels(site))
names(myColors) = unique(site)

names(myColors) = levels(site)
plot(scores[,1:2],pch=16,col=myColors[site],asp=1)
ordispider(scores[,1:2],group=site,col=myColors)
ordiellipse(scores[,1:2],group=site,col=myColors,draw="polygon",label=T)

# --------- constrained ordination (by site)

pp1=capscale(ibs~site)
axes2plot=c(1,2)
scores=data.frame(pp1$CCA$wa[,axes2plot])
ggplot() + geom_point(data=scores,aes(scores[,1],scores[,2],color=site),alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])

# ordination after taking AWAY effect of site (imagine site was a nuisance covariate)
scores=data.frame(pp1$CA$u[,axes2plot]) # after taking out effect of site
ggplot(scores,aes(scores[,1],scores[,2],color=site)) + geom_point(alpha=0.5)+theme_bw()+xlab(names(scores)[1])+ylab(names(scores)[2])


