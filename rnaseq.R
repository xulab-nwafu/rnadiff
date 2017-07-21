library(methods)
args <- commandArgs(trailingOnly = TRUE)
target_file <- args[1]
rnaseq_file <- args[2]
control_tag <- args[3]
case_tag    <- args[4]

# target_file <- "demo/targets.txt"
# rnaseq_file <- "demo/rnaseq.tsv"
# control_tag <- "wt"
# case_tag    <- "mu"

library(edgeRun)

# load data
targets <- readTargets(file = target_file)
x<-read.delim(file = rnaseq_file, row.names=1, stringsAsFactors=FALSE, as.is = TRUE, sep = "\t")

# Preprocess
y <- DGEList(counts=x, group=targets$Treatment)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# QC plot
pdf("edgeRun.dir/QC.pdf")
pairs(x, pch=".", upper.panel = NULL)
plotBCV(y)
plotMDS(y, method = "bcv")
dev.off()

# DEG
z <- UCexactTest(y, pair=c(control_tag, case_tag))
de <- decideTestsDGE(z, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(de)

# DE tags
detags <- rownames(y)[as.logical(de)]
uptags<-rownames(y)[as.logical(de>0)]
downtags<-rownames(y)[as.logical(de<0)]
write.table(detags, "edgeRun.dir/de.txt", quote=F,col.names = F, row.names=F)
write.table(downtags, "edgeRun.dir/de.down.txt", quote=F,col.names = F, row.names=F)
write.table(uptags, "edgeRun.dir/de.up.txt", quote=F,col.names = F, row.names=F)

# MA plot
pdf("edgeRun.dir/MA.pdf")
plotSmear(z, de.tags=detags)
dev.off()

# Diff table
diff <- topTags(z, n = nrow(y))
# Cpm table
expr <-cpm(y)[rownames(diff),]
# Summary table
sum<-cbind(diff$table, expr)
write.table(sum, "edgeRun.dir/rnaseq.sum.tsv", quote=F, col.names=NA, row.names = TRUE, sep="\t")

# Save image
save.image(file = "edgeRun.dir/rnaseq.RData")
