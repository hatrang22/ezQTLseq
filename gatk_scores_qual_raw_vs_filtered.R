library('gridExtra')
library('ggplot2')


args <- commandArgs(trailingOnly = TRUE)
if (length(args) !=4) {
  stop("At least 4 arguments must be supplies:\nvcf count table raw\nvcf count table filtered\nTitle for graphs (SNPs ...)\npdf output file\n", call.=FALSE)
}

#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

table1  = args[1] # vcf count table raw SNPs|indels
table2  = args[2] # vcf count table filtered SNPs|indels
title     = args[3] # Title for graphs ("SNPs" ...)
outfile   = args[4] # pdf output file

data1 <- read.csv(table1, header = T, na.strings=c("","NA"), sep = "\t")
data2 <- read.csv(table2, header = T, na.strings=c("","NA"), sep = "\t")

dim(data1)
dim(data2)
leg1 <- gsub(" ", "_", paste("Raw",title))
leg2 <- gsub(" ", "_", paste("Filtered",title))
VCF <- rbind(data1,data2)
VCF$Variant <- factor(c(rep(leg1, dim(data1)[1]),rep(leg2, dim(data2)[1])))

craw  <- '#A9E2E4'
cfilt <- '#F4CCCA'

QUAL<- ggplot(VCF, aes(x=QUAL, fill=Variant)) + geom_density(alpha=0.3) + scale_x_continuous(trans = 'log10')
DP <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3)
QD <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) + geom_vline(xintercept=20, size=0.7)
FS <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3)
MQ <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) + geom_vline(xintercept=40, size=0.7)
MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3)
SOR <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) #+ geom_vline(xintercept=c(4, 10), size=1, colour = c(craw,cfilt))
ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) #+
         #geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(craw,craw,cfilt,cfilt)) + xlim(-30, 30)


VCF <- rbind(data2)
VCF$Variant <- factor(c(rep(rep(leg2, dim(data2)[1]))))

QUALf<- ggplot(VCF, aes(x=QUAL, fill=Variant)) + geom_density(alpha=0.3) + scale_x_continuous(trans = 'log10')
DPf <- ggplot(VCF, aes(x=DP, fill=Variant)) + geom_density(alpha=0.3)
QDf <- ggplot(VCF, aes(x=QD, fill=Variant)) + geom_density(alpha=.3) + geom_vline(xintercept=20, size=0.7)
FSf <- ggplot(VCF, aes(x=FS, fill=Variant)) + geom_density(alpha=.3)
MQf <- ggplot(VCF, aes(x=MQ, fill=Variant)) + geom_density(alpha=.3) + geom_vline(xintercept=40, size=0.7)
MQRankSumf <- ggplot(VCF, aes(x=MQRankSum, fill=Variant)) + geom_density(alpha=.3)
SORf <- ggplot(VCF, aes(x=SOR, fill=Variant)) + geom_density(alpha=.3) #+ geom_vline(xintercept=c(4, 10), size=1, colour = c(craw,cfilt))
ReadPosRankSumf <- ggplot(VCF, aes(x=ReadPosRankSum, fill=Variant)) + geom_density(alpha=.3) #+
         #geom_vline(xintercept=c(-10,10,-20,20), size=1, colour = c(craw,craw,cfilt,cfilt)) + xlim(-30, 30)


pdf(outfile, height=40, width=15)
theme_set(theme_gray(base_size = 16))
grid.arrange(QUAL,QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum,QUALf,QDf, DPf, FSf, MQf, MQRankSumf, SORf, ReadPosRankSumf, nrow=9)
dev.off()


