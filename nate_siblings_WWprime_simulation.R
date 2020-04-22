sibs <- function(n, k, id) { #genotypes of full sibs
  freq <- rbeta(n, 1, 1)
  sire <- rbinom(n, size=2, prob=freq)
  dam <- rbinom(n, size=2, prob=freq)
  kids <- lapply(1:k, function(x) rhyper(n, sire, 2-sire, 1) + rhyper(n, dam, 2-dam, 1))
  geno <- Reduce(rbind, kids)
  rownames(geno) <- paste0("clutch", id, "_kid", 1:k)
  colnames(geno) <- paste0("locus", 1:ncol(geno))
  geno
}

clutch_1 <- sibs(10000, 10, 1)
clutch_2 <- sibs(10000,  5, 2)

W <- rbind(clutch_1, clutch_2)
#remove SNPs with no variance in sample
W <- W[,!apply(W, 2, function(x) length(unique(x))==1)]
A <- scale(W) %*% t(scale(W)) / ncol(W)
A[1:5,1:5]
library(ggplot2); library(reshape2)
ggplot(melt(A)) + geom_tile(aes(x=Var1, y=Var2, fill=value)) +
  scale_fill_gradient2() + theme(axis.text.x=element_text(angle=90)) +
  xlab("Individual") + ylab("Individual") 
?melt
