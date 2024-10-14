##可视化
#conda activate seurat
rm(list=ls())

library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
#library(ComplexHeatmap)
library(data.table)
library(pheatmap)
library(stringi)
library(dplyr)



## 提取 out_SCENIC.loom 信息
outDir <- '/Users/ceci/Project/LifeSpan/seurat/sample/S4varGenes3000_PCA25_Res0.6_K15/subset/MC/M1_varGenes3000_PCA10_Res0.4/pySCENIC/MC_All/'
scenicLoomPath='MC_All_out_SCENIC.loom'
loom <- open_loom(paste0(outDir,scenicLoomPath))
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")
n=t(scale(t( getAUC(regulonAUC[,] ))))
dim(n)
dim(regulonAUC)
regulonAUC[1:4,1:2]
sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
dim(sub_regulonAUC)

DF <- getAUC(sub_regulonAUC)
write.table(DF,file=paste0(outDir,'regulonActivity_by_cells.txt'),sep='\t',quote=F)
#python 生成celltypes分类


########biogroup########
## 提取 细胞分类信息 - 年龄组别的TF活性
biogroup <- read.table(paste0(outDir,'cells_by_BioGroup.txt'),sep='\t',header=T,row.names= 1)
head(biogroup)
selectedResolution <- "BioGroup"
cellsPerGroup <- split(rownames(biogroup), biogroup[,selectedResolution])
regulonActivity_by_biogroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
head(regulonActivity_by_biogroup)
write.table(regulonActivity_by_biogroup,file=paste0(outDir,'regulonActivity_by_BioGroup_Scaled.txt'),sep='\t',quote=F)


########biogroup########
regulonActivity_by_biogroup_Scaled <- t(scale(t(regulonActivity_by_biogroup),
                                          center = T, scale=T))
rss=regulonActivity_by_biogroup_Scaled
df = do.call(rbind,
            lapply(1:ncol(rss), function(i){
              dat= data.frame(
                path  = rownames(rss),
                cluster =   colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, median)
              )
            }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster)
n = rss[top5$path,]
pdf(file=paste0(outDir,'regulonActivity_by_BioGroup_Scaled_Top5.pdf'),width=3,height=7)
pheatmap(n,show_rownames = T, cellwidth = 8, cellheight = 8, fontsize = 8,cluster_cols=FALSE)
dev.off()

write.table(n,file=paste0(outDir,'regulonActivity_by_BioGroup_Scaled_Top5.txt'),sep='\t',quote=F)




#——————————————————————————————————s1-提取 细胞分类信息 - 各组别TF活性————————————————————————————————————#
########celltype########
cellTypes <- read.table(paste0(outDir,'cells_by_seurat_clusters.txt'),sep='\t',header=T,row.names= 1)
head(cellTypes)
selectedResolution <- "seurat_clusters"
cellsPerGroup <- split(rownames(cellTypes), cellTypes[,selectedResolution])
regulonActivity_by_cellTypes <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
head(regulonActivity_by_cellTypes)
write.table(regulonActivity_by_cellTypes,file=paste0(outDir,'regulonActivity_by_seurat_clusters_Scaled.txt'),sep='\t',quote=F)

#————————————————————————————————————————s2-各组别排名前几的TF画图——————————————————————————————————————————#
########celltype########
#regulonActivity_by_biogroup / regulonActivity_by_cellTypes / regulonActivity_by_biogroup_celltype
regulonActivity_by_cellTypes_Scaled <- t(scale(t(regulonActivity_by_cellTypes),
                                          center = T, scale=T))
rss=regulonActivity_by_cellTypes_Scaled
df = do.call(rbind,
            lapply(1:ncol(rss), function(i){
              dat= data.frame(
                path  = rownames(rss),
                cluster =   colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, median)
              )
            }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster)
n = rss[top5$path,]
pdf(file=paste0(outDir,'regulonActivity_by_seurat_clusters_Scaled_Top5.pdf'),width=3,height=6)
pheatmap(n,show_rownames = T, cellwidth = 8, cellheight = 8, fontsize = 8)
dev.off()
write.table(n,file=paste0(outDir,'regulonActivity_by_seurat_clusters_Scaled_Top5.txt'),sep='\t',quote=F)





########biogroup_celltype########
## 提取 细胞分类信息 - 年龄组别的TF活性
biogroup_celltype <- read.table(paste0(outDir,'cells_by_Seuratclusters_BioGroup.txt'),sep='\t',header=T,row.names= 1)
head(biogroup_celltype)
selectedResolution <- "Seuratclusters_BioGroup"
cellsPerGroup <- split(rownames(biogroup_celltype), biogroup_celltype[,selectedResolution])
regulonActivity_by_biogroup_celltype <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))
head(regulonActivity_by_biogroup_celltype)
write.table(regulonActivity_by_biogroup_celltype,file=paste0(outDir,'regulonActivity_by_Seuratclusters_BioGroup_Scaled.txt'),sep='\t',quote=F)


########biogroup_celltype########
regulonActivity_by_biogroup_celltype_Scaled <- t(scale(t(regulonActivity_by_biogroup_celltype),
                                          center = T, scale=T))
rss=regulonActivity_by_biogroup_celltype_Scaled
df = do.call(rbind,
            lapply(1:ncol(rss), function(i){
              dat= data.frame(
                path  = rownames(rss),
                cluster =   colnames(rss)[i],
                sd.1 = rss[,i],
                sd.2 = apply(rss[,-i], 1, median)
              )
            }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster)
n = rss[top5$path,]
pdf(file=paste0(outDir,'regulonActivity_by_Seuratclusters_BioGroup_Scaled_Top5.pdf'),width=15,height=25)
pheatmap(n,show_rownames = T, cellwidth = 5, cellheight = 5, fontsize = 5)
dev.off()
write.table(n,file=paste0(outDir,'regulonActivity_by_Seuratclusters_BioGroup_Scaled_Top5.txt'),sep='\t',quote=F)
