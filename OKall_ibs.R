# installing required packages: do this once, then remark with #

setwd('~/Dropbox/mega2019/mega2019_clean/RAD') # change this to where your scp'd files are
bams=read.table("bams")[,1] # list of bam files
bams=sub(".bam","",bams)
length(bams)

ma = as.matrix(read.table("OKall.ibsMat"))
dimnames(ma)=list(bams,bams)

ma[1:6,1:6]

hc=hclust(as.dist(ma),"ave")

plot(hc,cex=0.8) # clustering of samples by IBS (great to detect clones or closely related individuals)

#-------- pruning clonal replicates

cutclones=0.2
abline(h=cutclones,col="red")
grps=cutree(hc,h= cutclones)

# retaining a single representative of each clonal group
pruned=c();i=1
for (i in unique(grps)) {
	pruned[i]=names(grps[grps==i])[1]
}
length(pruned)
ma1=ma[pruned,pruned]
hc=hclust(as.dist(ma1),"ave")
plot(hc,cex=0.8) 

# also removing both K4 and O5 = onbviously there is some mixup going on with sample names
pruned=pruned[-which(pruned %in% c("K4","O5"))]
length(pruned)
write.table(paste(pruned,".bam",sep=""),file="bams_noclones",quote=F, col.names=F, row.names=F)

# scp bams_noclones to TACC, rerun angsd
