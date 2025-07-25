setwd('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/')
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(paletteer)
library(ggsci)
library(msigdbr)
# Negation of %in%
'%!in%' <- Negate('%in%')
getPalette2 = colorRampPalette(pal_npg("nrc", alpha = 1)(5))
colors = getPalette2(12)
prismatic::color(pal_npg("nrc", alpha = 1)(5))
## Path
Data_Dir = '/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/data/'

## Data Prepare ----
exp_all = read.table(paste0(Data_Dir, '/LSC_counts.txt'),header = T)
exp_all[1:4,1:4]
colnames(exp_all)
exp_all$maxlen = exp_all$End-exp_all$Start## not match
save(exp_all,file = paste0(Data_Dir, '/Data_split_pre.Rdata'))

load(file = paste0(Data_Dir, '/Data_split_pre.Rdata'),verbose = T)
exp_count = as.data.frame(exp_all)
table(duplicated(exp_count$Geneid))
#FALSE  TRUE 
#36591    10 

exp_count <- exp_count[!duplicated(exp_count$Geneid),]
table(duplicated(exp_count$Geneid))
##FALSE 
##36591 
table(exp_count$Geneid == '-')
##FALSE  TRUE 
##36591    0

rownames(exp_count) <- exp_count$Geneid
colnames(exp_count)[7:58] = str_split(colnames(exp_count)[7:58],'_',simplify = T)[,1]
sample_group = data.frame(row.names = colnames(exp_count)[7:58],
                          group = substr(x = colnames(exp_count)[7:58], 1, regexpr("\\d", colnames(exp_count)[7:58])-1), ## sub group name
                          treat = substr(x = colnames(exp_count)[7:58], str_length(colnames(exp_count)[7:58]), str_length(colnames(exp_count)[7:58])) ## sub group name
)
sample_group$Size = ifelse(str_detect(rownames(sample_group),'B'),'Large','Small')
sample_group$group = ifelse(str_detect(sample_group$group,'AML'),'AML',ifelse(str_detect(sample_group$group,'U'),'U937','KG1')
)
sample_group$treat = ifelse(str_detect(sample_group$treat,"\\d"),'Non',ifelse(str_detect(sample_group$treat,'J'),'JAS','LTA')
)
sample_group$sample_pair = gsub('[A-Z]$','',rownames(sample_group))
sample_group$ComTreat = paste0(sample_group$Size,'_',sample_group$treat)
unique(sample_group$ComTreat)
sample_group$sample_pairN = str_sub(sample_group$sample_pair,str_length(sample_group$sample_pair), str_length(sample_group$sample_pair))
sample_group$Group_ComTreat = paste0(sample_group$group,'_',sample_group$ComTreat)

exp_count <- exp_count[,rownames(sample_group)] # rm gene length

save(exp_count, sample_group,file = paste0(Data_Dir, '/Data_DE_pre.Rdata'))

## All sample PCA -----
load(file = paste0(Data_Dir, '/Data_DE_pre.Rdata'),verbose = T)
sample_needs = 'All_sample'
outDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/DEseq2_',sample_needs,'/')
if (!dir.exists(outDir)){
  dir.create(outDir)}

data1 = exp_count
sample_group$group = factor(sample_group$group,levels = c('AML',"KG1","U937"))
sample_group$treat = factor(sample_group$treat,levels = c('Non',"JAS","LTA"))
sample_group$Size = factor(sample_group$Size,levels = c('Large',"Small"))
sample_group$ComTreat = factor(sample_group$ComTreat,levels = c("Large_Non","Large_LTA","Small_Non","Small_JAS"))
sample_group$sample_pair = as.factor(sample_group$sample_pair)
levels(sample_group$sample_pair)
sample_group$sample_pairN = as.factor(sample_group$sample_pairN)
sample_group$Group_ComTreat = as.factor(sample_group$Group_ComTreat)

table(str_detect(rownames(data1),'^AC'))
table(str_detect(rownames(data1),'^AL'))

dds <- DESeqDataSetFromMatrix(countData=data1,colData=sample_group, design=~sample_pairN+ComTreat)
table(rowSums(data1)>=1)
table(rowSums(data1)>=10)

smallestGroupSize <- 4
keep <- rowSums(counts(dds, normalized=FALSE) >= 10) >= smallestGroupSize
dds <- dds[keep,]
table(keep)
#FALSE  TRUE 
#15722 20869
#### PCA -----
vst <- vst(dds)

plotPCA(vst,intgroup="Group_ComTreat")+stat_ellipse(level = 0.95)
data <- plotPCA(vst,intgroup="Group_ComTreat",returnData = TRUE)
plotPCA(vst,intgroup="Group_ComTreat")+
  stat_ellipse(level = 0.95)+
  labs(title=" ") +
  geom_text_repel(aes(PC1, PC2, label = rownames(data))) +
  theme(aspect.ratio = 1,
        axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=15),
        axis.title.x  = element_text(size=15),
        axis.ticks=element_blank(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "white", colour = "white"),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text = element_text(size=15),
        plot.background=element_blank())+ 
  #xlab(paste0( " PC1(79%)" ))+
  #ylab(paste0( " PC2(10%)" ))+
  scale_color_manual(values=colors)
ggsave(paste0(outDir,'/PCA_',sample_needs,'.pdf'),width = 8,height = 7)

## Sample Group split -----
## pair sample design = ~ subjects + condition. subjects are samples 

load(file = paste0(Data_Dir, '/Data_DE_pre.Rdata'),verbose = T)

theme_lsc = theme(aspect.ratio = 1,
                  axis.text.y   = element_text(size=12),
                  axis.text.x   = element_text(size=12),
                  axis.title.y  = element_text(size=15),
                  axis.title.x  = element_text(size=15),
                  axis.ticks=element_blank(),
                  panel.grid = element_blank(),
                  legend.title = element_blank(),
                  legend.text = element_text(size=15))

DEseq_spGroup <- function(count = exp_count, coldata = sample_group,
                          col_group = c("Large_Non"='#3C5488',"Large_LTA"='#4DBBD5',"Small_Non"='#F39B7F',"Small_JAS"='#E64B35'),
                          sample_rm = NULL,SampleGroup = 'AML'
) {
  sample_needs = SampleGroup
  outDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/DEseq2_',sample_needs,'/')
  if (!dir.exists(outDir)){
    dir.create(outDir)}
  
  coldata = coldata[coldata$group==sample_needs,]
  data1 = count[,rownames(coldata)] 
  print(dim(data1))
  coldata$treat = factor(coldata$treat,levels = c('Non',"JAS","LTA"))
  coldata$Size = factor(coldata$Size,levels = c('Large',"Small"))
  coldata$ComTreat = factor(coldata$ComTreat,levels = c("Large_Non","Large_LTA","Small_Non","Small_JAS"))
  coldata$sample_pair = as.factor(coldata$sample_pair)
  levels(coldata$sample_pair)
  coldata$sample_pairN = as.factor(coldata$sample_pairN)
  dds <- DESeqDataSetFromMatrix(countData=data1,colData=coldata, design=~sample_pairN+ComTreat)
  print(table(rowSums(data1)>=10))
  
  smallestGroupSize <- 8
  keep <- rowSums(counts(dds, normalized=FALSE) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  print(table(keep))
  
  #### PCA -----
  vst <- vst(dds)
  
  plotPCA(vst,intgroup="ComTreat")+stat_ellipse(level = 0.95)
  data <- plotPCA(vst,intgroup="ComTreat",returnData = TRUE)
  plotPCA(vst,intgroup="ComTreat")+
    stat_ellipse(level = 0.95)+
    labs(title=" ") +
    geom_text_repel(aes(PC1, PC2, label = rownames(data))) +
    theme(aspect.ratio = 1,
          axis.text.y   = element_text(size=12),
          axis.text.x   = element_text(size=12),
          axis.title.y  = element_text(size=15),
          axis.title.x  = element_text(size=15),
          axis.ticks=element_blank(),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.key = element_rect(fill = "white", colour = "white"),
          panel.background = element_rect(color = 'black', fill = 'transparent'),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.text = element_text(size=15),
          plot.background=element_blank())+ 
    #xlab(paste0( " PC1(79%)" ))+
    #ylab(paste0( " PC2(10%)" ))+
    scale_color_manual(values=col_group)
  
  ggsave(paste0(outDir,'/PCA_',sample_needs,'.pdf'),width = 8,height = 7)
  
  ### DE -----
  dds <- DESeq(dds)
  dds
  res = results(dds)
  res
  head(results(dds, tidy=TRUE))
  summary(res) #summary of results
  ## Sort summary list by p-value
  res <- res[order(res$padj),]
  head(res)
  
  save(dds,coldata,vst,file = paste0(outDir,'/DEseq2_',sample_needs,'_result.Rdata'))
  
  cut_off_pvalue = 0.05
  cut_off_logFC = 0.585
  
  #### Large_LTA_vs_Large_Non ----
  res1 <- results(dds,contrast = c('ComTreat','Large_LTA','Large_Non'))
  res1
  print(summary(res1)) 
  res1 <- res1[order(res1$padj),]
  head(res1)
  normD = as.data.frame(counts(dds,normalized=TRUE))
  result_data <- merge(as.data.frame(res1),normD[,rownames(coldata)[which(coldata$ComTreat%in%c('Large_LTA','Large_Non'))]],by='row.names',sort=FALSE)
  names(result_data)[1] <- 'Gene'
  #result_data = as.data.frame(res1)
  result_data$change = ifelse(result_data$pvalue < cut_off_pvalue & abs(result_data$log2FoldChange) >= cut_off_logFC, 
                              ifelse(result_data$log2FoldChange> cut_off_logFC ,'Up','Down'), 'NotChange')
  cat(paste0(sample_needs,':'))
  cat('\n')
  print(table(result_data$change))
  result_data_sign = result_data[result_data$change!='NotChange',]
  openxlsx::write.xlsx(result_data, file=paste0(outDir,'Large_LTA_vs_Large_Non_gene_all.xlsx'))
  openxlsx::write.xlsx(result_data_sign, file=paste0(outDir,'Large_LTA_vs_Large_Non_gene_signif.xlsx'))
  
  ggplot(result_data, 
         aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
    geom_point(size=3,alpha = 0.8, shape = 19,stroke = 0) +
    scale_color_manual(values=c('Down' = '#505A7F', 'Up' = '#AD5053', 'Not changed' = "#d2dae2"))+
    geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8) +
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      min.segment.length = 0,
      color = "black",
      segment.color = "black", show.legend = FALSE )+
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=-5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      color = "black",
      segment.color = "black",
      min.segment.length = 0
    )+
    labs(x="log2FoldChange",y="-log10(pvalue)",title = 'Large_LTA_vs_Large_Non')+
    theme_bw(  
    )+
    theme_lsc
  ggsave(paste0(outDir, 'Large_LTA_vs_Large_Non_volcano.pdf'),width = 7,height = 6)
  
  
  
  #### Small_JAS_vs_Small_Non ----
  res1 <- results(dds,contrast = c('ComTreat','Small_JAS','Small_Non'))
  res1
  print(summary(res1)) 
  res1 <- res1[order(res1$padj),]
  head(res1)
  normD = as.data.frame(counts(dds,normalized=TRUE))
  result_data <- merge(as.data.frame(res1),normD[,rownames(coldata)[which(coldata$ComTreat%in%c('Small_JAS','Small_Non'))]],by='row.names',sort=FALSE)
  names(result_data)[1] <- 'Gene'
  #result_data = as.data.frame(res1)
  result_data$change = ifelse(result_data$pvalue < cut_off_pvalue & abs(result_data$log2FoldChange) >= cut_off_logFC, 
                              ifelse(result_data$log2FoldChange> cut_off_logFC ,'Up','Down'), 'NotChange')
  cat(paste0(sample_needs,':'))
  cat('\n')
  print(table(result_data$change))
  result_data_sign = result_data[result_data$change!='NotChange',]
  openxlsx::write.xlsx(result_data, file=paste0(outDir,'Small_JAS_vs_Small_Non_gene_all.xlsx'))
  openxlsx::write.xlsx(result_data_sign, file=paste0(outDir,'Small_JAS_vs_Small_Non_gene_signif.xlsx'))
  
  ggplot(result_data, 
         aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
    geom_point(size=3,alpha = 0.8, shape = 19,stroke = 0) +
    scale_color_manual(values=c('Down' = '#505A7F', 'Up' = '#AD5053', 'Not changed' = "#d2dae2"))+
    geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8) +
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      min.segment.length = 0,
      color = "black",
      segment.color = "black", show.legend = FALSE )+
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=-5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      color = "black",
      segment.color = "black",
      min.segment.length = 0
    )+
    labs(x="log2FoldChange",y="-log10(pvalue)",title = 'Small_JAS_vs_Small_Non')+
    theme_bw(  
    )+
    theme_lsc
  
  ggsave(paste0(outDir, 'Small_JAS_vs_Small_Non_volcano.pdf'),width = 7,height = 6)
  
  
  
  #### Small_Non_vs_Large_Non ----
  res1 <- results(dds,contrast = c('ComTreat','Small_Non','Large_Non'))
  res1
  print(summary(res1)) 
  res1 <- res1[order(res1$padj),]
  head(res1)
  normD = as.data.frame(counts(dds,normalized=TRUE))
  result_data <- merge(as.data.frame(res1),normD[,rownames(coldata)[which(coldata$ComTreat%in%c('Small_Non','Large_Non'))]],by='row.names',sort=FALSE)
  names(result_data)[1] <- 'Gene'
  #result_data = as.data.frame(res1)
  
  result_data$change = ifelse(result_data$pvalue < cut_off_pvalue & abs(result_data$log2FoldChange) >= cut_off_logFC, 
                              ifelse(result_data$log2FoldChange> cut_off_logFC ,'Up','Down'), 'NotChange')
  cat(paste0(sample_needs,':'))
  cat('\n')
  print(table(result_data$change))
  result_data_sign = result_data[result_data$change!='NotChange',]
  openxlsx::write.xlsx(result_data, file=paste0(outDir,'Small_Non_vs_Large_Non_gene_all.xlsx'))
  openxlsx::write.xlsx(result_data_sign, file=paste0(outDir,'Small_Non_vs_Large_Non_gene_signif.xlsx'))
  
  ggplot(result_data, 
         aes(x = log2FoldChange, y = -log10(pvalue), colour=change)) +
    geom_point(size=3,alpha = 0.8, shape = 19,stroke = 0) +
    scale_color_manual(values=c('Down' = '#505A7F', 'Up' = '#AD5053', 'Not changed' = "#d2dae2"))+
    geom_vline(xintercept=0,lty=4,col="gray",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="gray",lwd=0.8) +
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      min.segment.length = 0,
      color = "black",
      segment.color = "black", show.legend = FALSE )+
    geom_text_repel(
      data = top_n(result_data[result_data$pvalue<0.05,], n=-5,log2FoldChange),
      aes(log2FoldChange, -log10(pvalue),label = Gene),
      size = 4.5,
      color = "black",
      segment.color = "black",
      min.segment.length = 0
    )+
    labs(x="log2FoldChange",y="-log10(pvalue)",title = 'Small_Non_vs_Large_Non')+
    theme_bw(  
    )+
    theme_lsc
  ggsave(paste0(outDir, 'Small_Non_vs_Large_Non_volcano.pdf'),width = 7,height = 6)
  
  
  
}

DEseq_spGroup(count = exp_count, coldata = sample_group,
              col_group = c("Large_Non"='#3C5488',"Large_LTA"='#4DBBD5',"Small_Non"='#F39B7F',"Small_JAS"='#E64B35'),
              SampleGroup = 'AML'
)
#keep
#FALSE  TRUE 
#18768 17823 

DEseq_spGroup(count = exp_count, coldata = sample_group,
              col_group = c("Large_Non"='#3C5488',"Large_LTA"='#4DBBD5',"Small_Non"='#F39B7F',"Small_JAS"='#E64B35'),
              SampleGroup = 'KG1'
)
#keep
#FALSE  TRUE 
#21355 15236

DEseq_spGroup(count = exp_count, coldata = sample_group,
              col_group = c("Large_Non"='#3C5488',"Large_LTA"='#4DBBD5',"Small_Non"='#F39B7F',"Small_JAS"='#E64B35'),
              SampleGroup = 'U937'
)
#keep
#FALSE  TRUE 
#21273 15318 
