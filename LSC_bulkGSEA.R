setwd('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/')
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggpubr)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(paletteer)
library(ggsci)

# Negation of %in%
'%!in%' <- Negate('%in%')
getPalette2 = colorRampPalette(pal_npg("nrc", alpha = 1)(5))
colors = getPalette2(12)
prismatic::color(pal_npg("nrc", alpha = 1)(5))

source('/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/gsea_plot2_modify.R')

## Path
#"Large_Non","Large_LTA","Small_Non","Small_JAS"

load(file = '/Users/liangtingting/hp/lab_file/WMM/RNAseq_CRH/Data/GSEA/human_geneset.Rdata')

GSEA_lscF = function(SampleGroup = 'AML',Experim = 'Large_LTA',Control = 'Large_Non'){
  dataDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/DEseq2_', SampleGroup)
  outDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/GSEA_', SampleGroup)
  
  if (!dir.exists(outDir)){
    dir.create(outDir)}
  
  geneList1 = openxlsx::read.xlsx(paste0(dataDir,'/', Experim, '_vs_', Control,'_gene_all.xlsx'))
  geneList = geneList1[,c('Gene','log2FoldChange')]
  
  ### check duplicated
  table(duplicated(geneList$Gene))
  
  rownames(geneList) <- geneList$Gene
  geneList_1 = arrange(geneList, desc(log2FoldChange))
  geneList1 <- geneList_1
  head(geneList1)
  
  table(geneList$log2FoldChange==0)
  
  id <- geneList1$log2FoldChange
  names(id) <- rownames(geneList1)
  head(id)
  
  
  egmt2 <- GSEA(id, TERM2GENE= all_sets[,c('gs_name','gene_symbol')] , 
                minGSSize = 10,
                eps = 0, 
                pvalueCutoff = 0.99,
                verbose=FALSE)
  
  gsea_results3 <- egmt2@result 
  #13250 obs; 11 vars
  #1553 obs; 11 vars
  #7151 obs; 11 vars
  
  write.xlsx(gsea_results3,file = paste0(outDir, '/', Experim, 'vs', Control,'Gsea_all_humanpathyway.xlsx'))
  
  ### |NES|>1，NOM pvalue<0.05，FDR（padj）<0.25
  g2 <- gsea_results3[gsea_results3$pvalue<0.05 & gsea_results3$p.adjust<0.25 & abs(gsea_results3$NES)>1,]
  g2 <- g2[order(g2$NES,decreasing = T),]
  #96 obs; 11 vars
  #1947 obs; 11 vars
  write.xlsx(g2,file = paste0(outDir, '/', Experim, 'vs', Control,'Gsea_humanpathyway_NES1_p0.05_FDR0.25.xlsx'))
  
  ###### save Rdata result 
  save(egmt2,file = paste0(outDir, '/', Experim, 'vs', Control,'Gsea_result.Rdata'))
  
} 


GSEA_lscF(SampleGroup = 'AML',Experim = 'Large_LTA',Control = 'Large_Non')
GSEA_lscF(SampleGroup = 'AML',Experim = 'Small_JAS',Control = 'Small_Non')
GSEA_lscF(SampleGroup = 'AML',Experim = 'Small_Non',Control = 'Large_Non')

GSEA_lscF(SampleGroup = 'KG1',Experim = 'Large_LTA',Control = 'Large_Non')
GSEA_lscF(SampleGroup = 'KG1',Experim = 'Small_JAS',Control = 'Small_Non')
GSEA_lscF(SampleGroup = 'KG1',Experim = 'Small_Non',Control = 'Large_Non')

GSEA_lscF(SampleGroup = 'U937',Experim = 'Large_LTA',Control = 'Large_Non')
GSEA_lscF(SampleGroup = 'U937',Experim = 'Small_JAS',Control = 'Small_Non')
GSEA_lscF(SampleGroup = 'U937',Experim = 'Small_Non',Control = 'Large_Non')

SelectGSEAPlot = function(SampleGroup = 'AML',Experim = 'Large_LTA',Control = 'Large_Non'){
  dataDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/GSEA_', SampleGroup)
  outDir = paste0(dataDir,'/GSEAFigure_',Experim, 'vs', Control)
  if (!dir.exists(outDir)){
    dir.create(outDir)}
  
  load(file = paste0(dataDir, '/', Experim, 'vs', Control,'Gsea_result.Rdata'))
  
  dataSel = openxlsx::read.xlsx(paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/Figure/', Experim, 'vs', Control,'Gsea_humanpathyway_NES1_p0.05_FDR0.25.xlsx'),
                                sheet = 2
  )
  dataSel$Change = ifelse(dataSel$NES>=1,'Up','Down')
  table(dataSel$Change)
  
  g3 = dataSel$Description
  
  for (i in 1:length(g3)) {
    
    pdf(paste0(outDir,'/',g3[i],'.pdf'),width = 6, height = 6)
    p = gseaplot2_mod_tt(egmt2,geneSetID = g3[i],
                         title = NULL,
                         color = "green",
                         color_hmp = c("#39518F","#7D2125"),
                         base_size = 10,
                         rel_heights = c(1.5, 0.5, 1),
                         subplots = 1:3,
                         pvalue_table = F,
                         fontface = "plain", 
                         pCol = "black", 
                         pvalSize = 6,
                         pvalX = 0.9, pvalY = 0.8,
                         ES_geom = "line"
    )
    print(p)
    dev.off()
    
  }
}

SelectGSEAPlot(SampleGroup = 'AML',Experim = 'Large_LTA',Control = 'Large_Non')
SelectGSEAPlot(SampleGroup = 'AML',Experim = 'Small_JAS',Control = 'Small_Non')
SelectGSEAPlot(SampleGroup = 'AML',Experim = 'Small_Non',Control = 'Large_Non')


#SampleGroup = 'AML';Experim = 'Large_LTA';Control = 'Large_Non'
#SampleGroup = 'AML';Experim = 'Small_JAS';Control = 'Small_Non'
SampleGroup = 'AML';Experim = 'Small_Non';Control = 'Large_Non'

dataDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/GSEA_', SampleGroup)
outDir = paste0(dataDir,'/GSEAFigure_',Experim, 'vs', Control)
if (!dir.exists(outDir)){
  dir.create(outDir)}

load(file = paste0(dataDir, '/', Experim, 'vs', Control,'Gsea_result.Rdata'))

dataSel = openxlsx::read.xlsx(paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/Figure/', Experim, 'vs', Control,'Gsea_humanpathyway_NES1_p0.05_FDR0.25.xlsx'),
                              sheet = 2
)
dataN = openxlsx::read.xlsx('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/Figure/Need_Gene.xlsx')
dataSel = dataSel[dataSel$Description%in%dataN[,paste0(Experim, '_', Control,'_P')],]
dataSel$Change = ifelse(dataSel$NES>=1,'Up','Down')
table(dataSel$Change)
##### bar plot ========
dataSel.long <-dataSel
dataSel.long$Description <- factor(dataSel.long$Description, levels = rev(dataSel.long$Description))
#dataSel.long[which(dataSel.long$Change == 'Down'), c('NES')] <- dataSel.long[which(dataSel.long$Change == 'Down'), c('NES')] * -1 

ggplot(dataSel.long,aes(NES,Description,fill=Change))+
  geom_col(width=0.6)+
  scale_fill_manual(values =  c('Down' = '#034e61', 'Up' = '#a00000'))+
  coord_fixed(ratio =0.5) +
  #coord_fixed(ratio =1) +
  #coord_fixed(ratio = 0.4) +
  theme_bw() +
  theme(legend.title = element_blank(),
        #aspect.ratio = 2.5,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+
  labs(x = "NES", y = "", title = 'GSEA')#+coord_flip() 

ggsave(file = paste0(outDir,"/SelectGSEA_",Experim,"vs", Control,".pdf"),width = 6.5,height = 2.8) #width = 10,height = 10  width = 8,height = 3.5 width = 9,height = 6


## New Select Pathway ====
#SampleGroup = 'AML';Experim = 'Large_LTA';Control = 'Large_Non'
#SampleGroup = 'AML';Experim = 'Small_JAS';Control = 'Small_Non'
SampleGroup = 'AML';Experim = 'Small_Non';Control = 'Large_Non'

dataDir = paste0('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/GSEA_', SampleGroup,'/')
outDir = paste0(dataDir,'/GSEAFigure_',Experim, 'vs', Control)
if (!dir.exists(outDir)){
  dir.create(outDir)}

load(file = paste0(dataDir, '/', Experim, 'vs', Control,'Gsea_result.Rdata'))

dataSel = openxlsx::read.xlsx(paste0(dataDir, Experim, 'vs', Control,'Gsea_humanpathyway_NES1_p0.05_FDR0.25.xlsx'),
                              sheet = 1
)
dataN = openxlsx::read.xlsx('/Users/liangtingting/hp/lab_file/ZMM/RNAseq_LSC/Figure/Need_Gene_v2.xlsx')
dataSel = dataSel[dataSel$Description%in%dataN[,paste0(Experim, '_', Control,'_P')],]
dataSel$Change = ifelse(dataSel$NES>=1,'Up','Down')
table(dataSel$Change)
##### bar plot ========
dataSel.long <-dataSel
dataSel.long$Description <- factor(dataSel.long$Description, levels = rev(dataSel.long$Description))
#dataSel.long[which(dataSel.long$Change == 'Down'), c('NES')] <- dataSel.long[which(dataSel.long$Change == 'Down'), c('NES')] * -1 

ggplot(dataSel.long,aes(NES,Description,fill=Change))+
  geom_col(width=0.6)+
  scale_fill_manual(values =  c('Down' = '#034e61', 'Up' = '#a00000'))+
  coord_fixed(ratio =0.5) +
  #coord_fixed(ratio =1) +
  #coord_fixed(ratio = 0.4) +
  theme_bw() +
  theme(legend.title = element_blank(),
        #aspect.ratio = 2.5,
        legend.text = element_text(size=14),
        axis.text = element_text(color="black",size=10),
        axis.title.y = element_text(color = "black",size=14),
        axis.title.x = element_text(color = "black",size=14))+
  labs(x = "NES", y = "", title = 'GSEA')#+coord_flip() 

ggsave(file = paste0(outDir,"/SelectGSEA_",Experim,"vs", Control,"New.pdf"),width = 8,height = 2.5) #width = 10,height = 10  width = 8,height = 3.5 width = 6.5,height = 2.8


