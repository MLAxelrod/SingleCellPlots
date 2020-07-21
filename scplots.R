###Single Cell RNA seq Plot Genes by Cell Type###
### Maggie Axelrod ###

#Libraries
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)
#################

### Dataset info ###
#Melanoma and HNSCC datasets processed per https://doi.org/10.1038/s41467-019-11738-0
#Healthy PBMC data available from 10X Genomics; cell types annotated per Seurat tutorial

#Replace "Path_to_your_files" with the path to your files
#Run these lines of code to import each processed dataset
mel<-readRDS("Path_to_your_files/melanoma.rds")
hnscc<-readRDS("Path_to_your_files/hnscc.rds")
pbmc<-readRDS("Path_to_your_files/pbmc.rds")

#Run the code below to make custom plotting functions that you can use to plot your gene of interest
melplt<-function(gene, ylab){
  VlnPlot(mel, features = c(gene), ncol = 1, group.by = "cellType")+
    scale_x_discrete(limits=c("Malignant", "Macrophage", "T cell", "B cell", "NK", "Endothelial", "CAF"))+
    scale_fill_manual(limits=c("Malignant", "Macrophage", "T cell", "B cell","NK", "Endothelial", "CAF" ),
                      values=c("gray30", "red", "blue", "darkgreen", "firebrick", "darkorange4", "darkorchid4"))+
    NoLegend()+ ylab(ylab)+ labs(title="Melanoma")+xlab("")
}

hplt<-function(gene, ylab){
  VlnPlot(hnscc, features = c(gene), ncol = 1, group.by = "cellType")+
    scale_x_discrete(limits=c("Malignant", "Macrophage", "T cell", "B cell", "Dendritic", "Mast", "Endothelial", "Fibroblast"))+
    scale_fill_manual(limits=c("Malignant", "Macrophage", "T cell", "B cell", "Dendritic", "Mast", "Endothelial", "Fibroblast"),
                      values=c("gray30", "red", "blue", "darkgreen", "firebrick", "lightsalmon4", "darkorange4", "darkorchid4"))+
    NoLegend()+ ylab(ylab)+ labs(title="HNSCC")+xlab("")
}


pbmcplt<-function(gene, ylab){
  VlnPlot(pbmc, features = c(gene), ncol = 1, group.by = "cellType")+
    scale_x_discrete(limits=c("CD14+ Mono", "FCGR3A+ Mono", "CD8 T", "Naive CD4 T", "Memory CD4 T", "B", "NK", "DC", "Platelet"))+
    scale_fill_manual(limits=c("CD14+ Mono", "FCGR3A+ Mono", "CD8 T", "Naive CD4 T", "Memory CD4 T", "B", "NK", "DC", "Platelet"),
                      values=c("gray30", "red", "blue", "darkgreen", "firebrick", "lightsalmon4", "darkorange4", "darkorchid4", "slateblue4"))+
    NoLegend()+ ylab(ylab)+ labs(title="Healthy PBMC")+xlab("")
}


###For each function, you need to input your gene and a y axis title
#Put your gene in quotes as the first argument
#what I have as second argument will italicize gene name, so need to put your gene in italic()

#examples
melplt("ACE2", expression(paste(italic(ACE2), " Expression Level")))
pbmcplt("ACE2", expression(paste(italic(ACE2), " Expression Level")))
hplt("ACE2", expression(paste(italic(ACE2), " Expression Level")))

pbmcplt("CD3E", expression(paste(italic(CD3E), " Expression Level")))
melplt("CD3E", expression(paste(italic(CD3E), " Expression Level")))
hplt("CD3E" , expression(paste(italic(CD3E), " Expression Level")))

#Heatmap of selected genes
DoHeatmap(mel, features=c("ACE2", "CD3E", "CD8A", "CD19", "CD14"), group.by = "cellType", size=3)+
  theme(text=element_text(size=13))


#Find markers differentially expressed betweeen macs and all others, and make volcano plot
Idents(mel)<-"cellType"
Idents(hnscc)<-"cellType"
T.mel<-FindMarkers(mel, ident.1 = "T cell") 
T.hnscc<-FindMarkers(hnscc, ident.1 = "T cell")
EnhancedVolcano(T.mel, lab=rownames(T.mel), x= 'avg_logFC', y= 'p_val_adj',
                FCcutoff = 0.5, ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top', legendLabSize = 10, legendIconSize = 2,
                drawConnectors = FALSE, labSize = 4, axisLabSize =11,
                title="Melanoma Dataset", subtitle = "Enriched in T cells", gridlines.major = FALSE, gridlines.minor = FALSE, caption="")

EnhancedVolcano(T.hnscc, lab=rownames(T.hnscc), x= 'avg_logFC', y= 'p_val_adj',
                FCcutoff = 0.5, ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendPosition = 'top', legendLabSize = 10, legendIconSize = 2,
                drawConnectors = FALSE, labSize = 4, axisLabSize =11,
                title="HNSCC Dataset", subtitle = "Enriched in T cells", gridlines.major = FALSE, gridlines.minor = FALSE, caption="")