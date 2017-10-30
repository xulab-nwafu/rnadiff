library(methods)
args <- commandArgs(trailingOnly = TRUE)
target_file <- args[1]
rnaseq_file <- args[2]
control_tag <- args[3]
case_tag    <- args[4]
out_dir     <- args[5]

# target_file <- "demo/rnaseq.conf"
# rnaseq_file <- "demo/rnaseq.tsv"
# control_tag <- "wt"
# case_tag    <- "mu"
# out_dir <- "edgeRun.dir"

library(edgeRun)

# load data
targets <- readTargets(file = target_file)
x<-read.delim(file = rnaseq_file, row.names=1, stringsAsFactors=FALSE, as.is = TRUE, sep = "\t")

# Check data
if  (!all(targets$Label == colnames(x)))
  stop("RNA-seq configuration file is not consistent with the data file",
       "\nRNA-seq conf label :\n", paste(targets$Label, "\t"),
       "\nRNA-seq data header:\n", paste(colnames(x), "\t")
       )
cmp_tags <- c(control_tag, case_tag)
if (!all(cmp_tags %in% unique(targets$Treatment)))
  stop("At least one tag for comparation is not defined in ", target_file, ":\n",
       "All defined tags are:\n", paste(unique(targets$Treatment), "\t"),
       "\nInput tags are:\n", paste(cmp_tags, "\t")
       )

# Preprocess
y <- DGEList(counts=x, group=targets$Treatment)
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

# QC plot
pdf(paste(out_dir, "/QC.pdf", sep=""))
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
write.table(detags, paste(out_dir, "/de.txt", sep=""), quote=F,col.names = F, row.names=F)
write.table(downtags, paste(out_dir, "/de.down.txt", sep=""), quote=F,col.names = F, row.names=F)
write.table(uptags, paste(out_dir, "/de.up.txt", sep=""), quote=F,col.names = F, row.names=F)

# MA plot
pdf(paste(out_dir, "/MA.pdf", sep=""))
plotSmear(z, de.tags=detags)
dev.off()

# Diff table
diff <- topTags(z, n = nrow(y))
# Cpm table
expr <-cpm(y)[rownames(diff),]
# Summary table
sum<-cbind(diff$table, expr)
write.table(sum, paste(out_dir, "/rnaseq.sum.tsv", sep=""), quote=F, col.names=NA, row.names = TRUE, sep="\t")

# Save image
save.image(file = paste(out_dir, "/rnaseq.RData", sep=""))
