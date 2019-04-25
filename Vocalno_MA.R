rm(list=ls())
library(DESeq2)
library(gplots)
library(ggplot2)
library(rgl)
library(rglwidget)
library(genefilter)
library(plot3D)
library(FactoMineR)
library(ggrepel)


temp<-read.csv("demo_table_volcano_MA.csv",header = T)

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
