rm(list=ls())
library(edgeRun)

# load data
targets <- readTargets("targets.txt")
x<-read.delim("rnaseq.tsv", row.names=1, stringsAsFactors=FALSE, as.is = TRUE)

# Preprocess
y <- DGEList(counts=x, group=targets$Treatment)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# QC plot
pdf("QC.pdf")
pairs(x, pch=".", upper.panel = NULL)
plotBCV(y)
plotMDS(y, method = "bcv")
dev.off()

# DEG
z <- UCexactTest(y,pair=c("wt","mu"))
de <- decideTestsDGE(z, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(de)

# DE tags
detags <- rownames(y)[as.logical(de)]
uptags<-rownames(y)[as.logical(de>0)]
downtags<-rownames(y)[as.logical(de<0)]
write.table(detags, "de.txt", quote=F,col.names = F, row.names=F)
write.table(downtags, "de.down.txt", quote=F,col.names = F, row.names=F)
write.table(uptags, "de.up.txt", quote=F,col.names = F, row.names=F)

# MA plot
pdf("MA.pdf")
plotSmear(z, de.tags=detags)
dev.off()

# Diff table
diff <- topTags(z, n = nrow(y))
# Cpm table
expr <-cpm(y)[rownames(diff),]
# Summary table
sum<-cbind(diff$table, expr)
write.table(sum, "rnaseq.sum.tsv", quote=F, col.names=NA, row.names = TRUE, sep="\t")

# Save image
save.image(file = "rnaseq.RData")
