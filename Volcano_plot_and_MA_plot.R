rm(list=ls())
library(gplots)
library(ggplot2)
library(rgl)
library(rglwidget)
library(genefilter)
library(plot3D)
library(FactoMineR)
dyn.load('/Library/Java/JavaVirtualMachines/jdk-9.0.1.jdk/Contents/Home/lib/server/libjvm.dylib')
library(xlsx)
library(ggrepel)

##############here I skip the steps of reading in the expression matrix and create condition matrix

#############Directly from the DESeq object construct
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~batch+condition) ########here I consider the batch effect in the model
dds=DESeq(dds,minReplicatesForReplace = 50)
norm=counts(dds,normalized = T)
rld=rlog(dds, blind=TRUE)
res=results(dds,contrast = c("condition","M","MC"))
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(norm), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
dim(subset(resdata, padj<0.05))
dim(subset(resdata, padj<0.01))
dim(subset(resdata, padj<0.05 & abs(log2FoldChange)>1))
dim(subset(resdata, padj<0.05 & abs(log2FoldChange)>0.5849625))
resdata<-resdata[,-c(4,5)]
gene<-readLines("target_gene.txt")#########this is the interested genes that you would like to label on the plot

temp<-resdata[,1:5]
temp$padj[which(is.na(temp$padj))]<-1 #####put all "NA" into 1 for adjusted pvalue
temp[which(is.na(temp$log2FoldChange)),]$log2FoldChange<-0 #######put all "NA" into 0 for log2FC
temp$label<-""  #####this is to save the label of your interested genes
temp$Significant<-"" #######significance shown on the volcano plot
temp$Regulation<-"" ########down or up regulation shown on the MA plot

temp$Regulation <- ifelse(temp$log2FoldChange >= 1, "Upregulated (FC > 2)", temp$Regulation)
temp$Regulation <- ifelse(temp$log2FoldChange <= -1, "Downregulated (FC < 0.5)", temp$Regulation)
temp$Regulation <- ifelse(temp$log2FoldChange > -1 & temp$log2FoldChange < 1, "NS", temp$Regulation)

temp$Significant <- ifelse(temp$padj < 0.05 & (temp$log2FoldChange > 1 | temp$log2FoldChange < -1), "FDR < 0.05, |FC|>2", "Not Sig")
temp[which(temp[,1]%in%gene),]$label<-matrix(unlist(strsplit(temp[which(temp[,1]%in%gene),][,1],'[_]')),byrow=T,ncol=2)[,2]

pdf("Valcona_test.pdf",14,10)
p<-ggplot(temp,aes(x=log2FoldChange,y=-log10(padj),label=label))
p<-p+geom_point(aes(color = Significant), size = 0.8)+scale_color_manual(values = c("red", "black"))
p<-p+theme_bw(base_size = 16) +geom_label_repel(point.padding = unit(0.1, "lines"))+labs(x = "log2FoldChange(M/MC)", y = "-log10(adjusted p-value)", title = "M -> MC", color = "")
p
dev.off()

pdf("MAplot.pdf",16,9)
p<-ggplot(temp,aes(x=log2(baseMean+1),y=log2FoldChange,label=label))
p<-p+geom_point(aes(color = Regulation), size = 0.8)+scale_color_manual(values = c("blue", "grey","red"))
p<-p+geom_hline(yintercept = c(0, -log2(2), log2(2)), linetype = c(1, 2, 2), color = c("black", "black", "black"))
p<-p+theme_bw(base_size = 16)+geom_label_repel(force = 1,box.padding = unit(0.35, "lines"),point.padding = unit(0.3,"lines"))
p<-p+labs(x = "Log2 mean expression", title = "M -> MC", color = "")
p
dev.off()
