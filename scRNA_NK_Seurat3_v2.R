Args <- commandArgs()
print(Args)

library('Seurat')
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(clustree)
library(viridis)

varGenes <- as.integer(Args[6])
CCAdim <- as.integer(Args[7])
PCAdim <- as.integer(Args[8])
#Resolution:
ResN <- Args[9]


PATH <- './seuratinteg'

outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_Res',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/NKcombined')

ctrl.data <- Read10X(data.dir='./cellranger/Cgroup/')
drug.data <- Read10X(data.dir='./cellranger/Agroup/')

#control
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Large", min.cells = 5)
ctrl$stim <- "Large"
ctrl[['percent.mt']] <- PercentageFeatureSet(object=ctrl, pattern='^MT-')
pdf(file=paste0(outDir,'/QC_ctrl.pdf'))
VlnPlot(ctrl,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
dev.off()
ctrl <- subset(ctrl, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = varGenes)
top10 <- head(VariableFeatures(ctrl), 10)
pdf(file=paste0(outDir,'/VariableFeatures_ctrl.pdf'))
plot1 <- VariableFeaturePlot(ctrl)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()

#drug
drug <- CreateSeuratObject(counts = drug.data, project = "Small", min.cells = 5)
drug$stim <- "Small"
drug[['percent.mt']] <- PercentageFeatureSet(object=drug, pattern='^MT-')
pdf(file=paste0(outDir,'/QC_drug.pdf'))
VlnPlot(drug,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3)
dev.off()
drug <- subset(drug, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)
drug <- NormalizeData(drug, verbose = FALSE)
drug <- FindVariableFeatures(drug, selection.method = "vst", nfeatures = varGenes)
top10 <- head(VariableFeatures(drug), 10)
pdf(file=paste0(outDir,'/VariableFeatures_drug.pdf'))
plot1 <- VariableFeaturePlot(drug)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off()


#Integration
NK.anchors <- FindIntegrationAnchors(object.list = list(ctrl, drug), dims = 1:CCAdim)
NK.combined <- IntegrateData(anchorset = NK.anchors, dims = 1:CCAdim)

DefaultAssay(NK.combined) <- "integrated"

#Scale&PCA&Cluster
NK.combined <- ScaleData(NK.combined, verbose = FALSE)
NK.combined <- RunPCA(NK.combined, npcs = 40, verbose = FALSE)
print("Run PCA Done!")
pdf(file=paste0(outDir,'/ElbowPlot.PCA.pdf'))
ElbowPlot(object = NK.combined, ndims = 50)
dev.off()
NK.combined <- RunUMAP(NK.combined, reduction = "pca", dims = 1:PCAdim)
NK.combined <- FindNeighbors(NK.combined, reduction = "pca", dims = 1:PCAdim)

#resolution:
#saveRDS(NK.combined,paste0(outDir,'/NKcombined_no_Resolution.rds'))
obj <- FindClusters(NK.combined, resolution = seq(0.1,0.8,by=0.1))
pdf(file=paste0(outDir,'/Clustree_resolution.pdf'))
clustree(obj)
dev.off()

NK.combined <- FindClusters(NK.combined, resolution = as.numeric(ResN))

UMAP=as.data.frame(Embeddings(object = NK.combined, reduction = "umap"))
write.table(UMAP,file=paste0(FilePrefix,'.umap.txt'),sep='\t',quote=F)
print("Save UMAP Done!")


#Plot Dimplot
pdf(file=paste0(outDir,'/Dimplot_','condition.pdf'))
DimPlot(NK.combined, reduction = "umap", group.by = "stim")
dev.off()
pdf(file=paste0(outDir,'/Dimplot_','condition2.pdf'))
DimPlot(NK.combined, reduction = "umap", split.by = "stim")
dev.off()
pdf(file=paste0(outDir,'/Dimplot_','cluster.pdf'))
DimPlot(NK.combined, reduction = "umap", label = TRUE)
dev.off()


#Save MetaData
MetaData=as.data.frame(NK.combined@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)

#Marker Gene
pbmc.markers=FindAllMarkers(object= NK.combined,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)
print("marker genes")
head(pbmc.markers)
MarkerGene=as.data.frame(pbmc.markers %>% group_by(cluster) %>% top_n(50,avg_log2FC))
write.table(MarkerGene,file=paste(FilePrefix,'.30MarkerGene.wilcox.txt'),quote =F,sep='\t')


#SaveExpData
#Data=as.data.frame(as.matrix(GetAssayData(object = Aging.integrated)))
#write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

#VarGenes=VariableFeatures(object = Aging.integrated)
#VarGeneData=Data[VarGenes,]
#write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.Integrate.txt'),sep='\t',quote=F)


###save file
#NK.combined$Group <- paste(NK.combined$stim, NK.combined$seurat_clusters, sep = "_")
#Idents(NK.combined) <- "Group"
saveRDS(NK.combined,paste0(outDir,'/NKcombined.rds'))


#Marker Gene
#DefaultAssay(NK.combined) <- "RNA"
#nk.markers <- FindConservedMarkers(NK.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
#write.table(nk.markers,file=paste0(outDir,'/ClusterMarkers_0.txt'),quote =F,sep='\t')


#DEGs
#Idents(NK.combined) <- "Group"
#Group <- FindMarkers(NK.combined, ident.1 = "DRUG_0", ident.2 = "CTRL_0", verbose = FALSE)
#write.table(Group, file=paste0(outDir,'/DEG_0.txt'),quote =F,sep='\t')


pdf(file=paste0(Dir,'VlnPlot_CD38.pdf'))
VlnPlot(combined.obj.AML, features = 'CD38', pt.size=0)+NoLegend()+theme(text = element_text(size = 8))
dev.off()
gene <- 'CD38'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = combined.obj.AML, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()


####################
###LSC6 score###
####################
library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(grid)
library(RColorBrewer)
library("UpSetR")
library(forcats)
library(Seurat)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(data.table)
library(grid)
library(RColorBrewer)

LSC_signature <- list(c("DNMT3B","CD34","GPR56","SOCS2","SPINK2","FAM30A"))
#FAM30A not in combined.obj.AML@assays$RNA@data
filtered_data <- combined.obj.AML@assays$RNA@data[which(rownames(combined.obj.AML@assays$RNA@data)%in%unlist(LSC_signature)),]
rownames(filtered_data)
Elsayed_model <- function(x){
  x[1]*0.0171+x[2]*0.109+x[3]*0.141+x[4]*0.0516+x[5]*0.189}
cells_scores <- apply(filtered_data,2,function(x)Elsayed_model(x))

Dir <- '/Users/xixi/Project/ZMM_scRNA/seuratinteg_3500/varGenes2000_CCA15_PCA15_Res0.5/LSC/'
#VlnPlot
pdf(file=paste0(Dir,'VlnPlot_Elsayed_LSC6_score.pdf'),width=5, height=2)
VlnPlot(All,features = "Elsayed_LSC6_score",pt.size=0) +NoLegend()
dev.off()

#FeaturePlot
pdf(file=paste0(Dir,'FeaturePlot_Elsayed_LSC6_score.pdf'),width=4.5, height=4)
FeaturePlot(combined.obj.AML,features = "Elsayed_LSC_score", pt.size = 0.6)+scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()

#ggplot
df1 <- data.frame(clusters=Idents(All),Elsayed_LSC6_score=All$Elsayed_LSC6_score)
df2 <- df1 %>% group_by(clusters) %>% summarise(median=median(Elsayed_LSC6_score))
quantiles <- quantile(df2$median,probs = seq(0, 1, 0.1))

pdf(file=paste0(Dir,'BoxPlot_Elsayed_LSC6_score.pdf'),width=4, height=3)
ggplot(df1,aes(x=fct_reorder(clusters, Elsayed_LSC6_score, .fun = median, .desc = TRUE),y=Elsayed_LSC6_score, fill=clusters)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  NoLegend()  +
coord_cartesian(ylim=c(-0.1,0.7)) +
  geom_hline(yintercept = quantiles[10], linetype="dashed")
dev.off()

####################
###LSC48 score###
####################
LSC48_score <- list(c("NPAL2",
"RBPMS",
"TRAF3IP2",
"PPP1R10",
"ATP1B1",
"NF1",
"RBPMS",
"FLJ13197",
"RBPMS",
"ABCG1",
"CLN5",
"LRRC8B",
"FRMD4B",
"ZFP30",
"C17orf86",
"C16orf5",
"TGIF2",
"RABGAP1",
"PPIG",
"GPR56",
"EIF2S3",
"NAB1",
"LRRC61",
"ATP1B1",
"ZNF500",
"CSDE1",
"C2CD2",
"PAQR6",
"FAM119B",
"ARPP-19",
"SETDB1",
"ZBTB39",
"RBPMS",
"SLC9A7",
"MAP3K7",
"ARL3",
"ZNF304",
"LOC552889",
"VGLL4",
"UBR5",
"PTCD2",
"CRKRS",
"IQGAP2",
"PLCH1",
"ARFGEF1",
"MAP3K7",
"PNPLA4"))
All <- AddModuleScore(object = All,features = LSC48_score,name = 'LSC48_score')

df1 <- data.frame(clusters=Idents(All),LSC48_score=All$LSC48_score1)
df2 <- df1 %>% group_by(clusters) %>% summarise(median=median(LSC48_score))
quantiles <- quantile(df2$median,probs = seq(0, 1, 0.1))

pdf(file=paste0(Dir,'BoxPlot_LSC48_score.pdf'),width=4, height=3)
ggplot(df1,aes(x=fct_reorder(clusters, LSC48_score, .fun = median, .desc = TRUE),y=LSC48_score, fill=clusters)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  NoLegend()  +
coord_cartesian(ylim=c(-0.2,0.25)) +
  geom_hline(yintercept = quantiles[10], linetype="dashed")
dev.off()


####################
#Predict the class of the cells using the markers and the expression of the BM cells form Van_Galen paper
####################
load(file="/Users/xixi/Project/ZMM_scRNA/DownloadData/VanGalen.obj.Rdata")
dt.list <- unlist(list(c(combined.obj.AML,merge.object)))
anchors <- FindTransferAnchors(reference = merge.object, query = combined.obj.AML, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = merge.object$predictionRF, dims = 1:30)
prediction <- factor(predictions$predicted.id,levels=c("HSC","Prog","GMP","ProMono","Mono","cDC","pDC","earlyEry","lateEry","ProB","B","Plasma","T","CTL","NK"))
names(prediction) <- rownames(predictions)
combined.obj.AML <- AddMetaData(object=combined.obj.AML, metadata = prediction, col.name = "prediction")
table(combined.obj.AML$orig.ident,combined.obj.AML$prediction)
aux_df <- data.frame(Condtion=combined.obj.AML$orig.ident,Predicted_cell_type=combined.obj.AML$prediction,Clusters=Idents(combined.obj.AML))
#Plot the celltypes by group
pdf(file=paste0(Dir,'ggplot_VanGalen_group.pdf'))
ggplot(aux_df, aes(Predicted_cell_type,fill=Condtion)) +
  geom_bar()
dev.off()
#Plot the celltypes by cluster
pdf(file=paste0(Dir,'ggplot_VanGalen_cluster.pdf'))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colors = getPalette(15)
ggplot(aux_df, aes(Clusters,fill=Predicted_cell_type)) +
  geom_bar() +
  scale_fill_manual(values = colors)
dev.off()
#Print the UMAPs with the prediction from BM
pdf(file=paste0(Dir,'UMAP_cluster.pdf'))
DimPlot(combined.obj.AML, reduction = "umap",  pt.size = 0.8) #+ scale_color_manual(values = colors)
dev.off()
pdf(file=paste0(Dir,'UMAP_VanGalen_cluster.pdf'))
DimPlot(combined.obj.AML, reduction = "umap",  pt.size = 0.8, label = TRUE, group.by = "prediction") +
  scale_color_manual(values = colors)
dev.off()


########AddModuleScore########
All <- AddModuleScore(object = All,features = GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS,name = 'GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS')
gene <- 'GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS1'
pdf(file=paste0(Dir,gene,'.pdf'))
FeaturePlot(object = All, features = gene)+
scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdYlBu")))
dev.off()
i<-"GO_REACTIVE_OXYGEN_SPECIES_METABOLIC_PROCESS1"
pdf(file=paste0(Dir,i,".Group.pdf"),width=5, height=2.5)
VlnPlot(All, features = c(i), slot = "data",  pt.size = 0, cols=c("#A379A2", "#D4D4D4","#00749F", "#CE5545","#F6AF15", "#81B5A1"))+NoLegend()+theme(text = element_text(size = 12))
dev.off()


#----------------------------------------------------------------------------
# cell cycle analysis
#----------------------------------------------------------------------------
cc.genes  <- readLines(con = "/Users/xixi/Project/ZMM_scRNA/DownloadData/scAML-main/tables/regev_lab_cell_cycle_genes.txt")
s.genes   <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
# Assign Cell-Cycle Scores
lsc_copy <- CellCycleScoring(object = lsc_copy, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

DefaultAssay(All) <- "integrated"
pdf(file=paste0(Dir,'DotPlot_CD7.pdf'),width=4, height=5)
DotPlot(All, features = 'CD7', cols = c("#A379A2", "#D4D4D4","#00749F", "#CE5545","#F6AF15", "#81B5A1"), dot.scale = 8) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()


#----------------------------------------------------------------------------
# 更改seurat的idents排列顺序
#----------------------------------------------------------------------------
levels(x=All）
levels(x=All) <- c('Blast','DC-like','Ery','GMP','HSPC','LSC')

#----------------------------------------------------------------------------
# 各亚群marker gene dotplot
#----------------------------------------------------------------------------
markers.to.plot <- c('LYZ','S100A9','VCAN','S100A8','CCL2',
'LTB','CCR7','HLA-DQB1','LSP1','DUSP5',
'HBD','HBA2','IGLL1','PRTN3','SLC40A1',
'CTSG','ELANE','HMGB3','MPO','FABP5',
'CD34','HLA-DQA2','HLA-DPB1','HLA-DPA1','IFITM3','CD74',
'LGALS3BP','SPINK2','CD7','ALDH1A1','ICAM3')

DefaultAssay(All) <- "RNA"
pdf(file=paste0(Dir,'DotPlot_Marker_RNA.pdf'),width=10.5, height=3.8)
DotPlot(All, features = rev(markers.to.plot), cols = c("#0078AA", "#D61C4E"), dot.scale = 8) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

DefaultAssay(All) <- "integrated"
pdf(file=paste0(Dir,'DotPlot_Marker_integrated.pdf'),width=10.5, height=3.8)
DotPlot(All, features = rev(markers.to.plot), cols = c("#0078AA", "#D61C4E"), dot.scale = 8) +
geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
