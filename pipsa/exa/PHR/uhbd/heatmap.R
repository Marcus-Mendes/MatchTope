library(gplots)
library(cluster)
library(pvclust)
pdf(file="matrix_2_h.pdf",width=9,height=9,paper="a4")
n=2
x=as.matrix(read.table("matrix"))
write.csv(x, file='table.csv')
x.dist=dist(x,method="maximum")
title="Electrostatic Distance maximum"
x.clust=hclust(x.dist,method="complete")
x.dend=as.dendrogram(x.clust)
x.clust.members=cutree(x.clust,k=n)
x.clust.coltypes=rainbow(n,s=0.4,start=0/6,end=5/6)
x.clust.colors=1:ncol(x)
order=order.dendrogram(x.dend)
for(count in 1:ncol(x)) {
   x.clust.colors[count]=x.clust.coltypes[x.clust.members[count]]
}
##### Snippet taken from SLmisc package to calculate constant colors for data
col=rainbow(512,s=0.8,v=1,start=0/6,end=5/6)
coln = col
margin=c(5,5)
#####
heatmap.2(as.matrix(x.dist),margins=c(14,14),col=coln,scale="none",distfun=function(n) n,Colv=x.dend,Rowv=x.dend,hclustfun=function(x) hclust(x.dist,method="complete"),trace="none",revC=TRUE,symkey=FALSE,denscol="white",density.info='density',densadj=1,ColSideColors=x.clust.colors,main="Electrostatic Distance maximum")
rm(x)
results <- pvclust(read.table("matrix"), method.hclust='complete', method.dist='cor', nboot=10000, parallel=FALSE)
plot(results)
output <- pvrect(results, alpha=0.97)
dev.off()
