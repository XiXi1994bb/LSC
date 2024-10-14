library(monocle)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggplot2)

All <- readRDS('/Users/xixi/Project/ZMM_scRNA/seuratinteg_3500/varGenes2000_CCA15_PCA15_Res0.5/NKcombined.rds')

NewPbmc = subset(All, subset=seurat_clusters!=11)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=6)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=13)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=8)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=12)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=9)
NewPbmc = subset(NewPbmc, subset=seurat_clusters!=0)

Idents(NewPbmc) <- "seurat_clusters"
NewPbmc<- RenameIdents(NewPbmc, `2` = "GMP", `4` = "GMP", `5` = "GMP", `1` = "Blast", `3` = "Blast", `7` = "Blast",`10` = "LSC")
NewPbmc$celltype <- Idents(NewPbmc)


DefaultAssay(NewPbmc) <- "RNA"
exprs=as.matrix(GetAssayData(object =NewPbmc[["RNA"]], slot = "counts"))
fd=as.data.frame(row.names(exprs))
colnames(fd)[1] <- "gene_short_name"
rownames(fd)<-row.names(exprs)

HSMM=newCellDataSet(exprs,featureData = new("AnnotatedDataFrame", data = fd),expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1) #
print(head(fData(HSMM)))
expressed_genes <- row.names(subset(fData(HSMM),num_cells_expressed >= 10))
print(head(pData(HSMM)))

pData(HSMM)$Total_mRNAs <- Matrix::colSums(HSMM@assayData$exprs)
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))

HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > 2000 &pData(HSMM)$Total_mRNAs < 10000]
HSMM <- detectGenes(HSMM, min_expr =0.02)

NewPbmc=NormalizeData(object=NewPbmc,normalization.method = 'LogNormalize',scale.factor=10000)
NewPbmc=FindVariableFeatures(NewPbmc,nfeatures = 2000,mean.cutoff = c(0.05, 6),dispersion.cutoff = c(0.5, Inf))


gene <- read.table('/Users/xixi/Project/ZMM_scRNA/seuratinteg_3500/varGenes2000_CCA15_PCA15_Res0.5/monocle/KEGG_CELL_CYCLE.v7.5.1.txt',sep='\t',header=1)
YGene=gene$gene
ordering_genes=VariableFeatures(object = NewPbmc)
for (g in VariableFeatures(object = NewPbmc)){
    if (g %in% YGene){
        print(g)
        ordering_genes=ordering_genes[ordering_genes!=g]
    }
}
VariableFeatures(object = NewPbmc)=ordering_genes


HSMM=setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM<- reduceDimension(HSMM,num_dim =20,max_components = 2,method = 'DDRTree') #降维
HSMM <- orderCells(HSMM,reverse=T)

Data <- NewPbmc@meta.data
Data2 <- Data[which(row.names(Data) %in% row.names(pData(HSMM))),]
pData(HSMM)$stim <- Data2$stim
pData(HSMM)$seurat_clusters <- Data2$seurat_clusters
pData(HSMM)$celltype <- Data2$celltype


Dir <- '/Users/xixi/Project/ZMM_scRNA/seuratinteg_3500/varGenes2000_CCA15_PCA15_Res0.5/monocle/'

pdf(file=paste0(Dir,"Monocle_pseudotime_removeCC.pdf"))
print("Save Pseudotime.pdf")
plot_cell_trajectory(HSMM,color_by="Pseudotime",cell_size=0.6)+
scale_color_viridis_c()
dev.off()


colorTrajStateFinal <- c("#D4D4D4", "#00749F", "#CE5545")
colorTrajStateFinal <- colorTrajStateFinal[c(1,2,3)]
pdf(file=paste0(Dir,"Monocle_celltype_removeCC.pdf"))
plot_cell_trajectory(HSMM, color_by = "celltype", cell_size=0.6) + scale_color_manual(values = colorTrajStateFinal)
dev.off()

colorTrajStateFinal <- c("#8CC0DE", "#FFD9C0")
colorTrajStateFinal <- colorTrajStateFinal[c(1,2)]
pdf(file=paste0(Dir,"Monocle_stim_removeCC.pdf"))
plot_cell_trajectory(HSMM, color_by = "stim", cell_size=0.6) + scale_color_manual(values = colorTrajStateFinal)
dev.off()
