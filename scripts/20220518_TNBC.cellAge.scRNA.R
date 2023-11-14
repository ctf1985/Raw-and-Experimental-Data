
  dir.create("scripts")
  dir.create("PDFs")
  dir.create("files")
  dir.create("origin_datas/TCGA",recursive = T)
  dir.create("origin_datas/scRNA",recursive = T)
  dir.create("results/files",recursive = T)
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("results/pdf",recursive = T)
  dir.create("results/ana",recursive = T)
}
options(stringsAsFactors = F)
source('/pub1/data/mg_projects/projects/codes/mg_base.R')
#1、单细胞聚类降维和细胞定义#####
dir.create('results/ana/01.scRNA')
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)

filenames=list.files( "origin_datas/scRNA/GSE176078/",full.names = F)
datalist=list()
for (i in 1:length(filenames)){
  filename=filenames[i]
  print(filename)
  dir.10x = paste0("origin_datas/scRNA/GSE176078/",filename)
  my.data <- Read10X(data.dir = dir.10x,gene.column = 1) 
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = filename,min.cells = 3, min.features = 250)
  datalist[[i]]$Sample=stringr::str_split_fixed(filename,'_',2)[,2]
}
rm(my.data)
names(datalist)=filenames

for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}

#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])

#细胞数统计
raw_meta=sce@meta.data
raw_count <- table(raw_meta$Sample)
raw_count
sum(raw_count)#38582

Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',feature2 = 'nCount_RNA',group.by = 'Sample')+
  theme(legend.position = 'none')
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nCount_RNA',group.by = 'Sample')+
  theme(legend.position = 'none')
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',feature2 = 'nFeature_RNA',group.by = 'Sample')


Feature_ber<-mg_merge_plot(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))

pearplot_befor<-VlnPlot(sce,group.by ='Sample', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0, 
                        ncol = 4)
pearplot_befor

#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset = nFeature_RNA > 100 & 
              nFeature_RNA < 6000 & 
              nCount_RNA > 100 )
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_meta=sce@meta.data

clean_count <- table(clean_meta$Sample)
clean_count
sum(clean_count)#38007
pearplot_after <- VlnPlot(sce,group.by ='Sample', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size = 0, 
                          ncol = 4)
pearplot_after

save(datalist,file = 'origin_datas/scRNA/datalist.RData')
#降维
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", 
                            nfeatures = 2000,
                            mean.cutoff=c(0.0125,3),
                            dispersion.cutoff =c(1.5,Inf))
sce <- ScaleData(sce, features =  rownames(sce))
sce <- RunPCA(sce, features = VariableFeatures(sce)) 
DimPlot(sce, reduction = "pca") 
ElbowPlot(sce, ndims=50, reduction="pca") 

sce <- RunUMAP(sce, dims=1:35, reduction="pca")
pp=DimPlot(sce,group.by = 'Sample',
           reduction="umap", 
           label = "T", pt.size = 0.2,
           label.size = 4) +
  theme(axis.line = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank())
pp
ggsave(filename = 'results/ana/01.scRNA/raw_merge.umap.pdf',plot = pp,he=9,wi=9)

#去批次
load('origin_datas/scRNA/datalist.RData')
for (i in 1:length(datalist)){
  datalist[[i]]<-NormalizeData(datalist[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
  datalist[[i]]<-FindVariableFeatures(datalist[[i]], 
                                      selection.method = "vst", 
                                      nfeatures = 2000,
                                      mean.cutoff=c(0.0125,3),
                                      dispersion.cutoff =c(1.5,Inf))
}
datalist <- FindIntegrationAnchors(object.list = datalist, dims = 1:35,
                                   reduction = c("cca", "rpca")[1])
sce <- IntegrateData(anchorset = datalist, dims = 1:35)

#ScaleData
DefaultAssay(sce) <- "integrated"
sce <- ScaleData(sce, features = rownames(sce))

sce=FindVariableFeatures(sce, 
                         selection.method = "vst", 
                         nfeatures = 2000,
                         mean.cutoff=c(0.0125,3),
                         dispersion.cutoff =c(1.5,Inf))



#PCA降维，选择合适的拐点
sce <- RunPCA(sce, features = VariableFeatures(sce)) 
dimplot1 <- DimPlot(sce, reduction = "pca",group.by = 'Sample') 
elbowplot1 <- ElbowPlot(sce, ndims=50, reduction="pca") 
sc_pca <- dimplot1+elbowplot1
sc_pca


Dims <- 30
#Umap降维
sce <- RunUMAP(sce, dims=1:Dims, reduction="pca")
###tsne 降维
sce <- RunTSNE(sce, 
               dims=1:Dims, 
               reduction="pca",
               perplexity=30,
               max_iter=1000)


figs1=mg_merge_plot(Feature_ber,pearplot_befor,pearplot_after,sc_pca,
                    nrow = 4,labels = c('A','B','C','D'),ncol = 1,widths = c(1,1,1,2))
ggsave(filename = 'results/ana/01.scRNA/FigS1.pdf',plot = figs1,he=12,width = 10)

library(clustree)
sce <- FindNeighbors(sce, dims = 1:Dims)
sce <- FindClusters(
  object = sce,
  resolution = c(seq(.1,1,.1))
)
colnames(sce@meta.data)
clustree(sce@meta.data, prefix = "integrated_snn_res.")

pdf('results/ana/01.scRNA/clust.snn_res.pdf',he=15,wi=15)
clustree(sce@meta.data, prefix = "integrated_snn_res.")
dev.off()

colnames(sce@meta.data)

#聚类分析
Resolution <- 0.1
sce <- FindNeighbors(object = sce, dims = 1:Dims)
sce <- FindClusters(object = sce, resolution = Resolution)
#细胞注释
DefaultAssay(sce)='RNA'

#CD8 T :1
VlnPlot(sce,features = c('CD3D','CD8A','CD8B','GZMA'),pt.size = 0,group.by = 'seurat_clusters')
#CD4 T :6
VlnPlot(sce,features = c('CD3D','CD4','CCR6','CCR7'),pt.size = 0,group.by = 'seurat_clusters')
#Treg
VlnPlot(sce,features = c('CD3D','CD4','FOXP3','RORC'),pt.size = 0,group.by = 'seurat_clusters')
#NK
VlnPlot(sce,features = c('KLRF1','FGFBP2','KLRC1'),pt.size = 0,group.by = 'seurat_clusters')
#NKT
VlnPlot(sce,features = c('CD2','KLRF1','FGFBP2'),pt.size = 0,group.by = 'seurat_clusters')
#Macrophage:2,9
VlnPlot(sce,features = c('CD163','CD68','CD14'),pt.size = 0,group.by = 'seurat_clusters')
#Monocyte:3,10
VlnPlot(sce,features = c('S100A12','VCAN','FCN1','S100A8'),pt.size = 0,group.by = 'seurat_clusters')

#B cells: 5
VlnPlot(sce,features = c('CD19','CD79A','MS4A1'),pt.size = 0,group.by = 'seurat_clusters')
#Plasma cell:4
VlnPlot(sce,features = c('CD79A','JSRP1'),pt.size = 0,group.by = 'seurat_clusters')
#Mast cell
VlnPlot(sce,features = c('TPSAB1','CPA3'),pt.size = 0,group.by = 'seurat_clusters')
#epithelial cells (EPCAM),0
VlnPlot(sce,features = 'EPCAM',pt.size = 0,group.by = 'seurat_clusters')

#Fibroblast:3,8
VlnPlot(sce,features = c('ACTA2','FAP','PDGFRB','NOTCH3'),pt.size = 0,group.by = 'seurat_clusters')
#endothelial cells:7
VlnPlot(sce,features = c('PECAM1'),pt.size = 0,group.by = 'seurat_clusters')


feat_gene<-c('CD3D','CD8A','CD8B','GZMA','CD4','CD19','CD79A','MS4A1','JSRP1','EPCAM','ACTA2','FAP','PDGFRB','NOTCH3','S100A8','CD163','CD68','PECAM1')

cell_anno<-data.frame(seurat_clusters=c('0','1','2','3','4','5','6','7','8','9','10'),
                      cell_type=c('Epithelial cell','CD8 T','Macrophage','Fibroblast','Plasma cell','B cell','CD4 T','Endothelial cell','Fibroblast','Macrophage','Monocyte'))
sce$cell_type<-sce$seurat_clusters
for (i in 1:nrow(cell_anno)){
  sce$cell_type=gsub(paste0('^',cell_anno$seurat_clusters[i],'$'),as.character(cell_anno$cell_type[i]),sce$cell_type)
}
table(sce$cell_type)
length(feat_gene)

length(feat_gene)

pdf('results/ana/01.scRNA/FigS2.pdf',he=8,width =16)
FeaturePlot(sce,
            features = feat_gene,
            pt.size = 0.1,reduction = 'tsne',ncol = 6)
dev.off()
#cnv分析
save(sce,file = 'sce.RData')
# load('sce.RData')
# library(copykat)

# copykat.test <- copykat(rawmat=sce@assays$RNA@counts,
#                         id.type="S",
#                         cell.line="no",
#                         ngene.chr=5,
#                         #每个染色体中至少有 5 个基因来计算 DNA 拷贝数
#                         win.size=25,
#                         #每个片段至少取 25 个基因
#                         KS.cut=0.15,
#                         #0-1,值越大灵敏度越低
#                         sam.name="TNBC",
#                         #随意固定一个名称
#                         distance="euclidean",
#                         n.cores=30
#                         #并行计算
# )
# save(copykat.test,file = 'copykat.test.RData')
#marker基因的筛选
#寻找差异基因时的差异倍数
Logfc = 0.5
#差异基因时最小的表达比例
Minpct = 0.5
DefaultAssay(sce) <- "RNA"
Idents(sce)<-'cell_type'
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)
write.table(sce.markers,'results/files/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')
write.table(sce.markers,'results/ana/01.scRNA/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')

### 选择前5个marker基因
Top5 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)  

Top5 <- intersect(unique(Top5$gene),rownames(sce@assays$RNA@meta.features))

sc_marker_dotplot <- DotPlot(object = sce, features = Top5,cols=c("blue", "red"),scale = T)+ 
  RotatedAxis()+ ggtitle("Top 5 Marker Genes")+ 
  theme(plot.title = element_text(hjust = 0.5)) +xlab('')

sc_marker_dotplot

#绘图
library(randomcoloR)
allcolour <- distinctColorPalette(60) 
length(table(sce@active.ident))
mycolor = allcolour[1:length(table(sce$seurat_clusters))]

colnames(sce@meta.data)
fig1a = DimPlot(sce,cols=allcolour,group.by = 'Sample',
                reduction="tsne",
                label = "F", 
                pt.size = 0.2,
                label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')+guides(colour = guide_legend(ncol = 1))

fig1a

fig1b<-DimPlot(sce,cols=mycolor,group.by = 'seurat_clusters',
               reduction="tsne",
               label = "F", 
               pt.size = 0.2,
               label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
fig1b
fig1c = DimPlot(sce,cols=mycolor,group.by = 'cell_type',
                reduction="tsne",
                label = "F", 
                pt.size = 0.2,
                label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
fig1c

fig1abc<-mg_merge_plot(fig1a,fig1b,fig1c,
                        nrow = 1,ncol=3,
                        labels = c('A','B','C'),
                        widths = c(1.2,1,1.3))

summary_cells <- as.data.frame(cbind(raw_count,clean_count))
write.table(summary_cells,'results/ana/01.scRNA/cell_count.txt',quote = F,row.names = T,sep='\t')
table(sce$seurat_clusters)
cell_anno
cell_anno$cell_num=as.data.frame(table(sce$seurat_clusters))[,2]
write.table(cell_anno,'results/ana/01.scRNA/cell_anno.txt',quote = F,row.names = F,sep='\t')
#
Idents(sce)='seurat_clusters'
library("ggplot2")
sample_clust<-as.data.frame(prop.table(table(sce$cell_type,sce$Sample)))
colnames(sample_clust)<-c("cluster","Patient","proportion")
write.table(sample_clust,'results/ana/01.scRNA/sample_clust.txt',quote = F,row.names = T,sep='\t')

clust_freq<-as.data.frame(table(sce$cell_type))
colnames(clust_freq)=c('cluster','cell_num')
clust_freq=clust_freq[order(clust_freq$cell_num,decreasing = T),]
clust_freq$cluster=factor(clust_freq$cluster,levels = clust_freq$cluster)
sample_clust$cluster=factor(sample_clust$cluster,levels =clust_freq$cluster)

fig1e1<-ggplot(sample_clust,aes(x = cluster,y = proportion,fill=Patient))+
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  theme_bw() + 
  theme(axis.ticks.length = unit(0.1, 'cm'),
        legend.position = "left") +xlab('')+
  coord_flip()+scale_y_continuous(expand = expand_scale(mult = c(0, 0)))
fig1e1

fig1e2<-ggplot(clust_freq,aes(x = cluster,y = cell_num))+
  geom_bar(stat="identity",fill='green')+ggtitle("") +
  theme_bw() + 
  theme(axis.ticks.length = unit(0, 'cm'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +coord_flip()+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)))+ylim(0,max(clust_freq$cell_num)+10)
fig1e2

fig1e3<-ggpubr::ggarrange(fig1e1,fig1e2,nrow = 1,ncol = 2,widths = c(2,1))
#细胞衰老相关通路的富集得分
library(clusterProfiler)
cellage.pathway<-clusterProfiler::read.gmt('origin_datas/cellAge.pathway.gmt')
write.table(cellage.pathway,'results/ana/01.scRNA/cellage.pathway.txt',quote = F,row.names = F,sep='\t')
write.table(cellage.pathway,'results/files/pathway.txt',quote = F,row.names = F,sep='\t')

ssGSEAScore_by_genes=function(gene.exp,genes){
  #library('GSVA')
  #library(GSEABase)
  #all.list=list()
  gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0("101"),geneIds=unique(genes),GSEABase::SymbolIdentifier()) 
  
  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  #detach('package:GSVA')
  #detach('package:GSEABase')
  #row.names(ssGSEA.geneset)
  return(ssGSEA.geneset)
}
as.character(unique(cellage.pathway$ont))

# cellage.pathway.score<-data.frame()
# exp=as.matrix(sce@assays$RNA@counts)
# for (i in as.character(unique(cellage.pathway$ont))){
#   path1=cellage.pathway[cellage.pathway$ont==i,]
#   score=ssGSEAScore_by_genes(gene.exp = exp,genes=path1$gene)
#   rownames(score)=i
#   cellage.pathway.score=rbind.data.frame(cellage.pathway.score,score)
# }
# rm(exp)
# cellage.pathway.score=t(cellage.pathway.score)
#save(cellage.pathway.score,file = 'cellage.pathway.score.RData')
head(cellage.pathway.score)

#读取CNV的结果
copykat.test<-read.delim('TNBC_copykat_prediction.txt',sep='\t',header = T)
head(copykat.test)
table(copykat.test$copykat.pred)
rownames(copykat.test)=copykat.test$cell.names
copykat.test=copykat.test[rownames(sce@meta.data),]
#添加分组
sce <- AddMetaData(sce, copykat.test$copykat.pred,col.name = "copykat.pred")
sce$copykat.pred[is.na(sce$copykat.pred)]<-'Unknown'
table(sce$copykat.pred)
# aneuploid   diploid 
# 7709      30298  
fig1f<- DimPlot(sce,cols=mycolor,group.by = 'copykat.pred',
                reduction="tsne",
                label = "F", 
                pt.size = 0.2,
                label.size = 5) +
  theme(axis.line = element_line(size=0.1, colour = "black"), 
        #axis.text = element_blank(), 
        #axis.title = element_blank(),
        axis.ticks = element_blank()
  ) +ggtitle('')
#各个样本中cnv预测的肿瘤和正常的比例
sample_clust1<-as.data.frame(prop.table(table(sce$copykat.pred,sce$Sample)))
colnames(sample_clust1)<-c("copykat.pred","Patient","proportion")
write.table(sample_clust1,'results/ana/01.scRNA/copykat.pred_sample.txt',quote = F,row.names = T,sep='\t')

fig1g<-ggplot(sample_clust1,aes(x = Patient,y = proportion,fill=copykat.pred))+
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  theme_bw() + 
  theme(axis.ticks.length = unit(0.1, 'cm'),
        legend.position = "right") +xlab('')+
  scale_y_continuous(expand = expand_scale(mult = c(0, 0)))+theme(axis.text.x = element_text(angle = 30,hjust = 1))

fig1ef<-mg_merge_plot(fig1e3,fig1f,fig1g,labels = c('E','F','G'),nrow = 1,ncol = 3,widths = c(1.5,1,1))
fig1=mg_merge_plot(fig1abc,sc_marker_dotplot,fig1ef,nrow = 3,ncol = 1,heights = c(1,1.2,1),labels = c('','D',''))
ggsave(filename = 'results/ana/01.scRNA/Fig1.pdf',plot = fig1,he=15,wi=18)


cellage.pathway.score.group<-merge(data.frame(cell.names=rownames(sce@meta.data),
                                              cell_type=sce@meta.data$cell_type,
                                              copykat.pred=sce@meta.data$copykat.pred),
                                   data.frame(cell.names=rownames(cellage.pathway.score),
                                              cellage.pathway.score),
                                   by='cell.names')
rownames(cellage.pathway.score.group)=cellage.pathway.score.group$cell.names
cellage.pathway.score.group=cellage.pathway.score.group[,-1]
head(cellage.pathway.score.group)
cellage.pathway.score.group.mean=data.frame()
for (i in c('aneuploid','diploid')){
  cellage.pathway.score.group1=cellage.pathway.score.group[cellage.pathway.score.group$copykat.pred==i,-2]
  cellage.pathway.score.group.mean1=aggregate(.~cell_type,data=cellage.pathway.score.group1,mean)
  cellage.pathway.score.group.mean1$copykat.pred=i
  cellage.pathway.score.group.mean=rbind.data.frame(cellage.pathway.score.group.mean,cellage.pathway.score.group.mean1)
}
head(cellage.pathway.score.group.mean)
table(cellage.pathway.score.group.mean$copykat.pred)

library(pheatmap)
library(ComplexHeatmap)
library(circlize)

#Z-score
cellage.pathway.score.group.mean1=apply(cellage.pathway.score.group.mean[,2:18],2,function(x){return(mosaic::zscore(x))})
cellage.pathway.score.group.mean=data.frame(cell_type=cellage.pathway.score.group.mean$cell_type,
                                            copykat.pred=cellage.pathway.score.group.mean$copykat.pred,
                                            cellage.pathway.score.group.mean1)


#aneuploid
mat1=cellage.pathway.score.group.mean[cellage.pathway.score.group.mean$copykat.pred=='aneuploid',]
rownames(mat1)=mat1$cell_type
col_anno <-HeatmapAnnotation(copykat.pred= mat1$copykat.pred)
names(col_anno)='aneuploid'

mat1=as.matrix(mat1[,c(as.character(unique(cellage.pathway$ont)))])

mean_value=apply(mat1,1,function(x){return(mean(x))})
range(mean_value)
bar = HeatmapAnnotation(barplot = anno_barplot(x = as.numeric(mean_value),
                                               border = F, 
                                               axis = TRUE,
                                               gp = gpar(fill = "green",
                                                         col = "white"),
                                               bar_with = 0.8,
                                               baseline=0,ylim = c(-1.5,2)))

p1= Heatmap(matrix = t(mat1),
            name="GSVA score",
            cluster_columns = F,
            cluster_rows = F,top_annotation = bar,
            #left_annotation= sam_ann,
            col = colorRamp2(c(-2,0,2), c("blue", "white", "red")),
            show_column_names = T,
            show_row_names = F,show_heatmap_legend = F,
            bottom_annotation = col_anno, 
            show_row_dend = F)
p1
#diploid
mat2=cellage.pathway.score.group.mean[cellage.pathway.score.group.mean$copykat.pred=='diploid',]
rownames(mat2)=mat2$cell_type
col_anno <-HeatmapAnnotation(copykat.pred= mat2$copykat.pred)
names(col_anno)='diploid'

mat2=as.matrix(mat2[,c(as.character(unique(cellage.pathway$ont)))])
mean_value=apply(mat2,1,function(x){return(mean(x))})
range(mean_value)
bar = HeatmapAnnotation(barplot = anno_barplot(x = as.numeric(mean_value),
                                               border = F, 
                                               axis = TRUE,
                                               gp = gpar(fill = "green",
                                                         col = "white"),
                                               bar_with = 0.8,
                                               baseline=0,ylim = c(-1.5,2)))


p2= Heatmap(matrix = t(mat2),
            name="Row Z-score GSVA score",
            cluster_columns = F,
            cluster_rows = F,top_annotation = bar,
            col = colorRamp2(c(-2,0,2), c("blue", "white", "red")),
            show_column_names = T,
            show_row_names = T,show_heatmap_legend = T,
            bottom_annotation  = col_anno, 
            show_row_dend = T)
p2
fig2<-p1+p2
fig2
pdf('results/ana/01.scRNA/Fig2.pdf',height = 7,width = 12)
print(fig2)
dev.off()
save(sce,file = 'sce.RData')
rm(sce,my.data)
save.image('project_001.RData')

##第二部分#####################

options(stringsAsFactors = F)
source('/pub1/data/mg_projects/projects/codes/mg_base.R')
#2、Bulk-Seq数据的准备####
dir.create('results/ana/02.data.pre')
#2.1 METABRIC####
cbio_cli<-read.delim('origin_datas/METABRIC/brca_metabric/data_clinical_sample.txt',sep='\t',header = T,skip = 4)
table(cbio_cli$ER_STATUS)
table(cbio_cli$HER2_STATUS)

table(cbio_cli$PR_STATUS)
cbio_cli=cbio_cli[which(cbio_cli$ER_STATUS=='Negative' & cbio_cli$HER2_STATUS=='Negative' & cbio_cli$PR_STATUS=='Negative'),]
cbio_cli=cbio_cli[which(cbio_cli$ER_STATUS=='Negative' & cbio_cli$HER2_STATUS=='Negative' & cbio_cli$PR_STATUS=='Negative'),]

cbio_cli<-data.frame(Samples=cbio_cli$PATIENT_ID,
                     oncotree_code=cbio_cli$ONCOTREE_CODE,
                     Grade=cbio_cli$GRADE,
                     TMB=cbio_cli$TMB_NONSYNONYMOUS,
                     Stage=cbio_cli$TUMOR_STAGE)
cbio_cli1<-read.delim('origin_datas/METABRIC/brca_metabric/data_clinical_patient.txt',sep='\t',header = T,skip = 4)
cbio_cli1<-data.frame(Samples=cbio_cli1$PATIENT_ID,
                      OS.time=ceiling(cbio_cli1$OS_MONTHS*30),
                      OS=cbio_cli1$OS_STATUS,
                      Chemotherapy=cbio_cli1$CHEMOTHERAPY,
                      RFS=cbio_cli1$RFS_STATUS,
                      RFS.time=ceiling(cbio_cli1$RFS_MONTHS*30),
                      NPI=cbio_cli1$NPI,
                      cellularity=cbio_cli1$CELLULARITY,
                      hormone_therapy=cbio_cli1$HORMONE_THERAPY,
                      radio_therapy=cbio_cli1$RADIO_THERAPY,
                      claudin_subtype=cbio_cli1$CLAUDIN_SUBTYPE,
                      Age=cbio_cli1$AGE_AT_DIAGNOSIS)
cbio_cli=merge(cbio_cli,cbio_cli1,by='Samples')
table(cbio_cli$Grade)

cbio_cli$OS.time
cbio_cli=cbio_cli[cbio_cli$OS.time>0,]
table(cbio_cli$OS.time<=10*365)
cbio_cli=cbio_cli[cbio_cli$OS.time<=10*365,]

table(cbio_cli$OS)
cbio_cli$OS=ifelse(cbio_cli$OS=='0:LIVING',0,1)
table(cbio_cli$RFS)
cbio_cli$RFS=ifelse(cbio_cli$RFS=='0:Not Recurred',0,1)

cbio_dat<-read.delim('origin_datas/METABRIC/brca_metabric/data_mrna_agilent_microarray.txt',sep='\t',header = T,check.names = F)
cbio_dat[1:4,1:4]
cbio_dat=cbio_dat[,-2]

head(gene_symbol)
gene_symbol=gene_symbol[which(gene_symbol$TYPE=='protein_coding'),"SYMBOL"]
com_gene<-intersect(cbio_dat$Hugo_Symbol,gene_symbol)
com_sample=intersect(colnames(cbio_dat),cbio_cli$Samples)
cbio_dat=cbio_dat[cbio_dat$Hugo_Symbol %in% com_gene,c('Hugo_Symbol',com_sample)]
cbio_dat=crbind2DataFrame(cbio_dat)
cbio_dat=aggregate(.~Hugo_Symbol,cbio_dat,mean)
rownames(cbio_dat)=cbio_dat$Hugo_Symbol
cbio_dat=cbio_dat[,-1]
cbio_dat=cbio_dat[,com_sample]
dim(cbio_dat)
range(cbio_dat)
# 15674   171
rownames(cbio_cli)=cbio_cli$Samples
cbio_cli=cbio_cli[com_sample,]
write.table(cbio_cli,'results/ana/02.data.pre/cbio_cli.txt',quote = F,sep='\t',row.names = F)
write.table(cbio_dat,'results/ana/02.data.pre/cbio_dat.txt',quote = F,sep='\t',row.names = T)

boxplot(cbio_dat[,1:4])
cbio_exp=cbio_dat
#2.2 GSE58812####
GSE58812=getGEOExpData('GSE58812')
GSE58812_cli=GSE58812$Sample
GSE58812_cli=data.frame(Samples=GSE58812_cli$Acc,
                        OS=GSE58812_cli$death,
                        OS.time=GSE58812_cli$`os (days)`,
                        ER=GSE58812_cli$`er-ihc`,
                        Her2=GSE58812_cli$`her2-ihc`,
                        PR=GSE58812_cli$`pr-ihc`)

table(GSE58812_cli$ER,GSE58812_cli$Her2,GSE58812_cli$PR)
GSE58812_cli$OS
GSE58812_cli$OS.time

GSE58812_exp=GSE58812$Exp$GPL570_54675_Data_col1
GSE58812_exp[1:4,1:4]
gene_anno<-GSE58812$Anno$GPL570[,c(1,11)]
gene_anno=gene_anno[gene_anno$V11 !='',]
gene_anno[,3]=stringr::str_split_fixed(string =gene_anno$V11,pattern = ' /// ',n = 2)[,2]
gene_anno=gene_anno[which(gene_anno$V3 ==''),c(1,2)]
colnames(gene_anno)=c('probe','gene')
GSE58812_exp=merge(gene_anno,data.frame(probe=rownames(GSE58812_exp),GSE58812_exp),by='probe')
GSE58812_exp[1:4,1:4]
GSE58812_exp=GSE58812_exp[,-1]
GSE58812_exp=aggregate(.~gene,GSE58812_exp,mean)
rownames(GSE58812_exp)=GSE58812_exp$gene
GSE58812_exp=GSE58812_exp[,-1]
GSE58812_exp=GSE58812_exp[gene_symbol,]
GSE58812_exp=na.omit(GSE58812_exp)
GSE58812_exp=log2(GSE58812_exp+1)
dim(GSE58812_exp)
#16416   107
GSE58812_exp[1:4,1:4]
write.table(GSE58812_exp,'results/ana/02.data.pre/GSE58812_exp.txt',quote = F,row.names = T,sep='\t')
write.table(GSE58812_cli,'results/ana/02.data.pre/GSE58812_cli.txt',quote = F,row.names = F,sep='\t')

#tcga#####
tcga_cli<-read.delim('origin_datas/TCGA/Merge_4ef0fcf3c7bb21e1aa2b9a7a96c248be_clinical.txt',sep='\t',header = T)
tcga_cli=data.frame(Samples=tcga_cli$A0_Samples,
                    Age=tcga_cli$A17_Age,
                    Gender=tcga_cli$A18_Sex,
                    T.Stage=tcga_cli$A3_T,
                    N.Stage=tcga_cli$A4_N,
                    M.Stage=tcga_cli$A5_M,
                    Stage=tcga_cli$A6_Stage,
                    OS=tcga_cli$A8_New_Event,
                    OS.time=tcga_cli$A8_New_Event_Time,
                    ER=tcga_cli$breast_carcinoma_estrogen_receptor_status,
                    PR=tcga_cli$breast_carcinoma_progesterone_receptor_status,
                    Her2=tcga_cli$lab_proc_her2_neu_immunohistochemistry_receptor_status)

table(tcga_cli$ER)
table(tcga_cli$PR)
table(tcga_cli$Her2)
tcga_cli=tcga_cli[which(tcga_cli$ER=='Negative' & tcga_cli$PR=='Negative' & tcga_cli$Her2=='Negative'),]
tcga_cli=tcga_cli[which(tcga_cli$OS.time>0),]
dim(tcga_cli)
table(tcga_cli$OS)
tcga_cli$Samples=paste0(tcga_cli$Samples,'-01')
rownames(tcga_cli)=tcga_cli$Samples
table(tcga_cli$OS.time>=30)
table(tcga_cli$OS.time<=10*365)
#tcga_cli=tcga_cli[which(tcga_cli$OS.time>=30),]

tcga_tmp<-read.delim('origin_datas/TCGA/Merge_TCGA-BRCA_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
sample_T<-intersect(tcga_cli$Samples,colnames(tcga_tmp)[substr(colnames(tcga_tmp),14,15)=='01'])
sample_N<-colnames(tcga_tmp)[substr(colnames(tcga_tmp),14,15)=='11']
length(sample_T)
#113
length(sample_N)
#113
tcga_tmp_log2<-log2(tcga_tmp[,c(sample_T,sample_N)]+1)
dim(tcga_tmp_log2)
# 25483   226
tcga_cli=tcga_cli[sample_T,]
tcga.group=rbind.data.frame(data.frame(Samples=sample_T,Type='T'),
                            data.frame(Samples=sample_N,Type='N'))
write.table(tcga_tmp_log2,'results/ana/02.data.pre/tcga_tpm_log.txt',quote = F,row.names = T,sep='\t')
write.table(tcga_cli,'results/ana/02.data.pre/tcga_cli.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.group,'results/ana/02.data.pre/tcga.group.txt',quote = F,row.names = F,sep='\t')

#3、细胞衰老相关通路在肿瘤和正常组织中的GSEA分析#########
dir.create('results/ana/03.GSEA/TCGA',recursive = T)
dim(tcga_tmp_log2)
head(tcga.group)
writeMatrix(dat = tcga_tmp_log2[,tcga.group$Samples],outpath = 'results/ana/03.GSEA/tcga_tmp_log2.txt',row = T,header = T)
writeMatrix(dat = tcga.group,outpath = 'results/ana/03.GSEA/tcga.group.txt',row = F,header = T)

mg_RunGSEA_wtl<-function(mod=c('exp_group','exp_gene','rank')[1],exp_Path=NULL, sample_group_path=NULL, outFolder=NULL,gene=NULL, column=NULL,lower=50,upper=50, gmt_Path=c("KEGG",'GO_BP','GO_CC','GO_MF','reactome','HALLMARK','TF','other')[1],plot_svg=FALSE,top=10,min=5,max=5000,outLog=T,gmt_Path1){
  
  if(is.null(exp_Path)|is.null(mod)|is.null(outFolder)|is.null(gmt_Path)){
    return(NULL)
  }
  if(plot_svg){
    svg='true'
  }else{
    svg='false'
  }
  if(gmt_Path=='KEGG'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.kegg.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_BP'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.bp.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_CC'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.cc.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_MF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.mf.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='reactome'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.reactome.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='HALLMARK'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/h.all.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='TF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c3.tft.v7.0.symbols.gmt')
  }else if(gmt_Path=='other'){
    gmt_Path=gmt_Path1
  }
  if(file.exists(paste0(getwd(),'/',outFolder))){
    outFolder=paste0(getwd(),'/',outFolder)
  }else if(!file.exists(outFolder)){
    dir.create(outFolder)
    if(file.exists(paste0(getwd(),'/',outFolder))){
      outFolder=paste0(getwd(),'/',outFolder)
    }
  }
  
  if(file.exists(paste0(getwd(),'/',gmt_Path))){
    gmt_Path=paste0(getwd(),'/',gmt_Path)
  }
  if(file.exists(paste0(getwd(),'/',exp_Path))){
    exp_Path=paste0(getwd(),'/',exp_Path)
  }
  
  command=NULL
  if(mod=='exp_group'){
    if(!is.null(exp_Path)&!is.null(sample_group_path)&!is.null(outFolder)){
      if(file.exists(paste0(getwd(),'/',sample_group_path))){
        sample_group_path=paste0(getwd(),'/',sample_group_path)
      }
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_group '
                     ,exp_Path,' ',sample_group_path,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }else if(mod=='exp_gene'){
    if(!is.null(exp_Path)&!is.null(gene)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_gene '
                     ,exp_Path,' ',gene,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max,' ',lower,' ',upper) 
    }
  }else if(mod=='rank'){
    if(!is.null(exp_Path)&!is.null(column)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar rank '
                     ,exp_Path,' ',column-1,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }
  
  if(!is.null(command)){
    if(MG_Grobal_System=='win'){
      command=paste0(MG_Grobal_baseFolder,'/jre/bin/java -jar ',command)
    }else{
      command=paste0('java -jar ',command)
    }
    print(paste0('RunGSEA CMD:',command))
    logs=system(command, intern = !outLog, 
                ignore.stdout = FALSE, ignore.stderr = FALSE, 
                wait = TRUE, input = NULL, show.output.on.console = TRUE, 
                minimized = FALSE, invisible = TRUE)
    if(outLog){
      if(logs==0){
        print('Run GSEA succ')
      }else{
        print('Run GSEA error')
      }
    }else{
      print('Runed GSEA')
      print(logs)
      logs=logs[grep('######/',logs)]
      if(length(logs)==1){
        logs=unlist(strsplit(logs[1],'/'))
        if(length(logs)>1){
          return(logs[2:length(logs)])
        }
      }
    }
  }
  return(NULL)
}
#TCGA
mg_RunGSEA_wtl(mod = 'exp_group',exp_Path = 'results/ana/03.GSEA/tcga_tmp_log2.txt'
               ,sample_group_path = 'results/ana/03.GSEA/tcga.group.txt'
               ,outFolder = 'results/ana/03.GSEA/TCGA/'
               ,gmt_Path = 'other',outLog=F,gmt_Path1=paste0(getwd(),'/origin_datas/cellAge.pathway.gmt'))
#结果读取
tcga_GSEA=parseGSEAResult('results/ana/03.GSEA/TCGA/my_analysis.Gsea.1652939373189/')

tcga.sig<-tcga_GSEA$EnrichTable[tcga_GSEA$EnrichTable$NP<0.05,"Term"]
sig_pathway=tcga.sig
sig_pathway

fig3a<-plot_GSEA_By_nodes(tcga_GSEA,TermNames = tcga.sig)
#4、细胞衰老相关通路在肿瘤和正常组织中的ssGSEA分析####
ssGSEAScore_by_genes=function(gene.exp,genes){
  #library('GSVA')
  #library(GSEABase)
  #all.list=list()
  gs=GSEABase::GeneSet(setName='GeneSet', setIdentifier=paste0("101"),geneIds=unique(genes),GSEABase::SymbolIdentifier()) 
  
  gsc <- GSEABase::GeneSetCollection(list(gs))
  fl <- tempfile()
  GSEABase::toGmt(gsc, fl)
  cgeneset=GSEABase::getGmt(fl)
  ssGSEA.geneset <- GSVA::gsva(as.matrix(gene.exp), cgeneset,method='ssgsea',
                               min.sz=1, max.sz=5000, verbose=TRUE)
  #detach('package:GSVA')
  #detach('package:GSEABase')
  #row.names(ssGSEA.geneset)
  return(ssGSEA.geneset)
}
pathway.score<-function(exp,gene){
  pathway_score<-data.frame()
  for (i in unique(gene[,1])){
    print(i)
    gene_set=gene[gene[,1]==i,"gene"]
    score=ssGSEAScore_by_genes(exp,gene_set)
    rownames(score)=i
    pathway_score=rbind.data.frame(pathway_score,score)
  }
  return(t(pathway_score))
}

cellage.pathway<-clusterProfiler::read.gmt('origin_datas/cellAge.pathway.gmt')
as.character(unique(cellage.pathway$ont))
table(cellage.pathway$ont)
table(cellage.pathway$gene)
#TCGA数据集ssGSEA分析
tcga.cellage.score<-pathway.score(exp=tcga_tmp_log2,gene = cellage.pathway)
write.table(tcga.cellage.score,'results/ana/03.GSEA/tcga.cellage.score.txt',quote = F,row.names = T,sep='\t')
write.table(tcga.cellage.score,'results/files/tcga.cellage.score.txt',quote = F,row.names = T,sep='\t')
#差异分析
diff_gene1<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  dat1=t(dat[dat$cluster==gr[1],-1])
  dat2=t(dat[dat$cluster==gr[2],-1])
  #dat3=t(dat[dat$cluster==gr[3],-1])
  #pathway=unique(c(rownames(dat1),rownames(dat2),rownames(dat3)))
  pathway=unique(c(rownames(dat1),rownames(dat2)))
  p_vale=data.frame()
  for (i in pathway){
    # d1=c(as.numeric(dat1[i,]),
    #      as.numeric(dat2[i,])
    #      #,as.numeric(dat3[i,])
    #      )
    # d2=factor(c(rep('A',ncol(dat1)),
    #             rep('B',ncol(dat2))
    #             #,rep('C',ncol(dat3)))
    #             ))
    dd1=wilcox.test(as.numeric(dat1[i,]),as.numeric(dat2[i,]))$p.value
    #logfc=log2(mean(dat1[i,])/mean(dat2[i,]))
    p_vale=rbind.data.frame(p_vale,data.frame(gene=i,p.value=dd1
                                              #,logfc=logfc
    ))
  }
  return(p_vale)
}

tcga.pathway.p<-diff_gene1(dat=t(tcga.cellage.score[tcga.group$Samples,]),
                           group=tcga.group$Type)
tcga.pathway.p$lab=ifelse(tcga.pathway.p$p.value<0.001,'***',ifelse(tcga.pathway.p$p.value<0.01,'**',ifelse(tcga.pathway.p$p.value<0.05,'*','')))

write.table(tcga.pathway.p,'results/ana/03.GSEA/tcga.pathway.p.txt',quote = F,row.names = F,sep='\t')


tcga.heat<-t(tcga.cellage.score)[tcga.pathway.p$gene,]
rownames(tcga.heat)=paste0(tcga.pathway.p$gene,tcga.pathway.p$lab)

#注释
tcga.anno<-data.frame(Type=tcga.group$Type)
rownames(tcga.anno)=tcga.group$Samples
tcga.anno=tcga.anno[order(tcga.anno$Type),,drop=F]

fig3b<-pheatmap::pheatmap(tcga.heat[,rownames(tcga.anno)],
                          scale = 'row',
                          show_colnames = F,annotation_col = tcga.anno,cluster_cols = F,cluster_rows =T,
                          color = colorRampPalette(c("blue", "white","red"))(100),
                          annotation_names_row = F,
                          breaks = unique(c(seq(-1, 1, length=100))))
fig3b <- ggplotify::as.ggplot(fig3b)
fig3<-mg_merge_plot(fig3a,fig3b,nrow = 2,ncol = 1,labels = c('A','B'))
ggsave('results/ana/03.GSEA/Fig3.pdf',fig3,height = 12,width = 15)
#5、衰老相关分子亚型的构建####
#Senescence Subtype
dir.create('results/ana/04.subtype')
sig_pathway=tcga.sig
sig_pathway
sig_gene=unique(cellage.pathway[cellage.pathway$ont %in% sig_pathway,"gene"])
length(sig_gene)#253
#
sig_gene=intersect(sig_gene,rownames(tcga_tmp_log2))
length(sig_gene)#186

#TCGA
tcga_cellage.exp<-tcga_tmp_log2[sig_gene,tcga_cli$Samples]
dim(tcga_cellage.exp)
write.table(tcga_cellage.exp,'results/ana/04.subtype/tcga_cellage.exp.txt',quote = F,row.names = T,sep='\t')
#186 113
#单因素cox分析
cox_batch<-function(dat,time,event){
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}
mg_risksocre.sig<-function(dat,os,os.time){
  #dat=tcga_cellage.exp_for
  #os=tcga_cli$OS
  #os.time=tcga_cli$OS.time
  tcga_dat_m<-cbind.data.frame(os=os,os.time=os.time,dat)
  tcga.cox <- cox_batch(t(tcga_dat_m[,3:ncol(tcga_dat_m)]),
                        time = tcga_dat_m$os.time,
                        event = tcga_dat_m$os)
  
  return(tcga.cox)
}
tcga.sig.cox<-mg_risksocre.sig(dat = t(tcga_cellage.exp[,tcga_cli$Samples]),
                               os = tcga_cli$OS,os.time = tcga_cli$OS.time)
table(tcga.sig.cox$p.value<0.05)
tcga.sig.cox.fit=tcga.sig.cox[which(tcga.sig.cox$p.value<0.05),]
tcga.sig.cox.fit=tcga.sig.cox.fit[order(tcga.sig.cox.fit$HR),]
bioForest=function(rt=null,col){
  #读取输入文件
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #输出图形
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
pdf('results/ana/04.subtype/Fig4a.pdf',height = 5,width = 7)
bioForest(rt = tcga.sig.cox.fit,col=c('blue','red'))
dev.off()
#一致性聚类
library(ConsensusClusterPlus)
clusterAlg_use=c('pam','hc','km','kmdist')[2]
distance_use=c('pearson','spearman','euclidean','canberra')[4]
#TCGA
df_exp=as.matrix(tcga_cellage.exp[rownames(tcga.sig.cox.fit),])
#df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
#df_exp=t(scale(t(df_exp)))
df_exp=t(scale(t(df_exp),scale = F))
dim(df_exp)
clust_subtype_TCGA = ConsensusClusterPlus(df_exp
                                          , maxK = 10, reps = 500, pItem = 0.8
                                          , pFeature =1
                                          , title = "TCGA_subtype"
                                          , clusterAlg = clusterAlg_use
                                          , distance = distance_use
                                          , innerLinkage = 'ward.D2'
                                          , plot = "pdf"
                                          , seed = 123456
                                          , writeTable = T)

k=3
tcga.subtype=data.frame(Cluster=clust_subtype_TCGA[[k]]$consensusClass)
rownames(tcga.subtype)<-colnames(df_exp)
tcga.subtype$Cluster=paste0('MC',tcga.subtype$Cluster)
tcga.subtype[which(tcga.subtype$Cluster=='MC1'),1]='TMC2'
tcga.subtype[which(tcga.subtype$Cluster=='MC2'),1]='TMC3'
tcga.subtype[which(tcga.subtype$Cluster=='MC3'),1]='TMC1'

tcga.subtype$Cluster=gsub('TMC','clust',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
tcga.subtype$Samples=rownames(tcga.subtype)

ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster)
            , add_text = ''
            #,labs = c('clust1','clust2')
)

write.table(data.frame(Samples=rownames(tcga.subtype),
                       time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster),
            'results/ana/04.subtype/tcga.subtype.txt',quote = F,row.names = F,sep='\t')
write.table(data.frame(Samples=rownames(tcga.subtype),
                       time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                       , event = tcga_cli[rownames(tcga.subtype),]$OS
                       , tcga.subtype$Cluster),
            'results/files//tcga.subtype.txt',quote = F,row.names = F,sep='\t')

fig4h<-ggplotKMCox(data.frame(time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                              , event = tcga_cli[rownames(tcga.subtype),]$OS
                              , tcga.subtype$Cluster)
                   , add_text = ''
                   ,palette = ggsci::pal_lancet()(9)[c(1,4,3)]
                   ,labs = c('clust1','clust2','clust3')
                   ,title = 'TCGA'
)
fig4h

ggsave('results/ana/04.subtype/fig4h.pdf',fig4h,height = 7,width = 7)


#关键基因表达的热图
library(ComplexHeatmap)
pdf('results/ana/04.subtype/Fig4b.pdf',height = 6,width = 6)
Heatmap(as.matrix(t(scale(t(tcga_tmp_log2[rownames(tcga.sig.cox.fit),tcga.subtype$Samples]))))
        , name = "Expr"
        , cluster_rows = F
        , cluster_row_slices = T
        , row_title_gp = gpar(fill = ggsci::pal_lancet()(9)[c(2,1,5)])
        , show_row_dend = F
        , column_split = tcga.subtype$Cluster
        , cluster_columns = F
        , cluster_column_slices=T
        , show_column_dend = F
        , show_column_names = F
        , col = circlize::colorRamp2(c(-4, 0, 4), c('#3B4992FF', 'white', '#EE0000FF'))
        , column_title_gp = gpar(fill = (ggsci::pal_lancet()(9))[c(1,4,3)])
        , border = TRUE
)
dev.off()
#snv 突变
library(maftools)

tcga_maf=read.maf('origin_datas/TCGA/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf',isTCGA = T)
tcga.maf1=subsetMaf(tcga_maf,tsb=intersect(tcga_maf@data$Tumor_Sample_Barcode,tcga.subtype1$Tumor_Sample_Barcode))

tcga.maf<-read.maf(tcga.maf1@data,
                    isTCGA=T,clinicalData = 'results/ana/04.subtype/anno.txt')

pdf('results/ana/04.subtype/Fig4c.pdf',width = 10,height = 5)
oncoplot(maf=tcga_maf,genes = rownames(tcga.sig.cox.fit),
         sortByAnnotation = T)
dev.off()
#cnv的表达
get_CNV_Preprocess=function(df_cnv){
  library(org.Hs.eg.db)
  library(clusterProfiler)
  df_cnv$`Gene Symbol`=gsub("\\..*","",df_cnv$`Gene Symbol`)
  rownames(df_cnv)=df_cnv$`Gene Symbol`
  
  library(TCGAutils)
  aliquot_id_to_submitter_id=UUIDtoBarcode(colnames(df_cnv)[-c(1:3)]
                                           ,from_type = 'aliquot_ids')
  colnames(df_cnv)[-c(1:3)]=aliquot_id_to_submitter_id[,2]
  colnames(df_cnv)=substr(colnames(df_cnv),1,15)
  df_cnv=df_cnv[,-which(duplicated(colnames(df_cnv)))]
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  
  ensg=df_cnv$`Gene Symbol`
  
  idmap=bitr(ensg,fromType="ENSEMBL",toType="SYMBOL",OrgDb="org.Hs.eg.db")
  idmap=dplyr::distinct(idmap,ENSEMBL,.keep_all=TRUE)## 一个基因匹配到多个geneid,随机取一个
  
  cnv.inds=which(!is.na(idmap$SYMBOL)) ## 去掉没有注释到gene symbol的行
  idmap=idmap[cnv.inds,]
  df_cnv=df_cnv[idmap$ENSEMBL,]
  df_cnv$`Gene Symbol`=idmap$SYMBOL
  
  df_cnv=dplyr::distinct(df_cnv,`Gene Symbol`,.keep_all=TRUE)
  dim(df_cnv)
  
  # rownames(df_cnv)=df_cnv$`Gene Symbol`
  df_cnv=df_cnv[,-c(2:3)]
  return(df_cnv)
}

cnv.all=read.delim('origin_datas/TCGA/Merge_GeneLevelCopyNumber.txt',sep = '\t'
                   ,stringsAsFactors = F,header = T,check.names = F)
cnv.all=get_CNV_Preprocess(cnv.all)
dim(cnv.all)
cnv.all[1:4,1:5]

table(substr(colnames(cnv.all)[-1],13,15))

names(table(substr(colnames(cnv.all)[-1],13,15)))

length(which(substr(colnames(cnv.all)[-1],13,15)=="-01"))

tcga.cnv.sp.selected=colnames(cnv.all)[-1][(which(substr(colnames(cnv.all)[-1],13,15)=="-01"))]

cnv.all=cnv.all[,c("Gene Symbol",tcga.cnv.sp.selected)]
dim(cnv.all)

table(substr(colnames(cnv.all)[-1],13,15))

get_CNV_Freq=function(df_cnv,genes_custom,genes_type=NULL){
  df_cnv=reshape2::melt(df_cnv,id.vars='Gene Symbol',measure.vars=colnames(cnv.all)[-1]
                        ,variable.name='Sample',value.name='CopyNum')
  head(df_cnv)
  cnv_frq=as.data.frame(matrix(data=0,nrow=1,ncol=5))
  colnames(cnv_frq)<-c('type','gene','CNV_gain','CNV_loss','none_CNV')
  r<-1
  for(i in 1:length(genes_custom))
  {
    cnv_gene<-df_cnv[which(df_cnv$`Gene Symbol` == genes_custom[i]),]
    total_counts<-dim(cnv_gene)[1]
    if(is.null(genes_type)){
      cnv_frq[r,1]<-NA
    }else{
      cnv_frq[r,1]<-genes_type[i]
    }
    cnv_frq[r,2]<-genes_custom[i]
    
    cnv_frq[r,3]<-(dim(cnv_gene[which(cnv_gene$CopyNum >0),])[1])/total_counts ## 在整个拷贝数数据集中
    cnv_frq[r,4]<-(dim(cnv_gene[which(cnv_gene$CopyNum <0),])[1])/total_counts ## CNV_loss
    cnv_frq[r,5]<-(total_counts/total_counts)-cnv_frq[r,3]-cnv_frq[r,4] ## none_CNV
    r<-r+1
  }
  cnv_frq<-cnv_frq[order(cnv_frq$CNV_loss),]
}

cnv.all[1:4,1:5]
table(substr(colnames(cnv.all[,-1]),14,15))
cnv.all=cnv.all[,intersect(c("Gene Symbol",tcga.subtype$Samples),colnames(cnv.all))]

cnv_freq=get_CNV_Freq(df_cnv=cnv.all,genes_custom = rownames(tcga.sig.cox.fit))

library(tidyr)
cnv_freq_1 <- cnv_freq %>%
  pivot_longer(c(CNV_gain, CNV_loss,none_CNV), names_to = "CNV", values_to = "value")

cnv_freq_1$gene = factor(cnv_freq_1$gene, levels=cnv_freq$gene)
fig4f=ggplot(cnv_freq_1,aes(x=gene,y=value,fill=CNV))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+scale_fill_manual(values = ggsci::pal_lancet('lanonc')(9)[c(2,1,3)])+
  labs(x="",y="Frequency",fill="CNV")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="top")
fig4f

get_Exp_Compare=function(df_exp,genes_custom=pathway.gene$gene){
  df_exp=data.frame(geneSymbol=rownames(df_exp),df_exp,stringsAsFactors = F,check.names = F)
  df_exp=reshape2::melt(df_exp,id.vars='geneSymbol',measure.vars=colnames(df_exp)[-1]
                        ,variable.name='Sample',value.name='Expr')
  head(df_exp)
  df_exp$sampleType=substr(df_exp$Sample,14,15)
  mode(df_exp$sampleType)='integer'
  
  sampleType.codes=read.table('/pub1/data/mg_projects/users/lishuang/public/TCGA/TCGA_SampleType_Codes.txt',header = T,sep='\t',stringsAsFactors = F,check.names = F)
  head(sampleType.codes)
  
  df_exp$cancer_type=NA
  for(code in unique(df_exp$sampleType)){
    inds=which(df_exp$sampleType==code)
    df_exp$cancer_type[inds]=sampleType.codes[which(sampleType.codes$Code==code),3]
  }
  head(df_exp)
  df_exp=df_exp[,c(1,2,3,5)]
  if(!is.null(genes_custom)){
    df_exp=df_exp[df_exp$geneSymbol %in% genes_custom,]
    rownames(df_exp)=NULL
  }
  return(df_exp)
}
dim(tcga_tmp_log2)
tcga.exp.cmp=get_Exp_Compare(tcga_tmp_log2,genes_custom = rownames(tcga.sig.cox.fit))
table(tcga.exp.cmp$cancer_type)
head(tcga.exp.cmp)
library(ggpubr)
tcga.exp.cmp$geneSymbol=factor(tcga.exp.cmp$geneSymbol,levels = com_gene)
fig4g<-ggboxplot(tcga.exp.cmp, x="geneSymbol", y="Expr", color = "cancer_type", 
                 palette = ggsci::pal_lancet('lanonc')(9)[c(1,2)]
                 # , facet.by = "fa"
                 , short.panel.labs = F)+
  stat_compare_means(aes(group=cancer_type), label = "p.signif", method = "t.test")+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))+xlab('')+
  ylab('log2(Gene Expression +1)')
fig4g
fig4fg<-mg_merge_plot(fig4f,fig4g,nrow = 1,ncol = 2,labels = c('F','G'))
ggsave('results/ana/04.subtype/Fig4fg.pdf',fig4fg,height = 5,width = 12)

#6、衰老相关亚型在相关通路中差异####
dir.create('results/ana/05.pathway.exp')
library(ComplexHeatmap)
library(circlize)

sig_pathway
head(cellage.pathway)
head(tcga_cellage.exp)
diff_pathway<-function(dat=tcga_cellage.exp[,tcga.subtype$Samples],group=tcga.subtype$Cluster){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  pathway=colnames(dat)[-1]
  p_vale=data.frame()
  for (i in pathway){
    dat1=dat[,c('cluster',i)]
    colnames(dat1)=c('cluster','gene')
    dd1=kruskal.test(dat1$gene,dat1$cluster)$p.value
    if(nrow(p_vale)==0){
      p_vale=data.frame(gene=i,p.value=dd1)
    }else{
      p_vale=rbind.data.frame(p_vale,data.frame(gene=i,p.value=dd1))
    }
  }
  return(p_vale)
}

tcga_cellage.exp.diff<-diff_pathway(dat = tcga_cellage.exp[,tcga.subtype$Samples],group=tcga.subtype$Cluster)
head(tcga_cellage.exp.diff)
tcga_cellage.exp.diff$lab<-ifelse(tcga_cellage.exp.diff$p.value<0.001,'***',ifelse(tcga_cellage.exp.diff$p.value<0.01,'**',ifelse(tcga_cellage.exp.diff$p.value<0.05,'*','')))
head(tcga_cellage.exp.diff)
tcga.pathway.exp<-tcga_cellage.exp[intersect(cellage.pathway[cellage.pathway$ont %in% sig_pathway,]$gene,
                                             rownames(tcga_cellage.exp)),]

tcga.pathway.exp_res=merge(data.frame(gene=tcga_cellage.exp.diff$gene,
                             lab=paste0(tcga_cellage.exp.diff$lab,tcga_cellage.exp.diff$gene)),
                  data.frame(gene=rownames(tcga.pathway.exp),tcga.pathway.exp),by='gene')
head(tcga.pathway.exp_res)
#rownames(tcga.pathway.exp_res)=tcga.pathway.exp_res$lab
#tcga.pathway.exp_res=tcga.pathway.exp_res[,-c(1,2)]
#head(tcga.pathway.exp_res)
colnames(tcga.pathway.exp_res)=gsub('\\.','-',colnames(tcga.pathway.exp_res))
library(pheatmap)
anno_col=data.frame(Cluster=tcga.subtype$Cluster)
rownames(anno_col)=tcga.subtype$Samples
anno_col=anno_col[order(anno_col$Cluster),,drop=F]
anno_row=cellage.pathway[which(cellage.pathway$ont %in% sig_pathway),]
anno_row=anno_row[anno_row$gene %in% tcga.pathway.exp_res$gene,]
ggsci::pal_lancet()(9)[c(1,4,3)]
#样本分组
sam_anno<-data.frame(Cluster=tcga.subtype$Cluster)
rownames(sam_anno)=tcga.subtype$Samples
sam_anno=sam_anno[order(sam_anno$Cluster),,drop=F]
sam_ann = rowAnnotation(df = sam_anno,
                        col = list(Cluster = c("clust1" ='#00468BFF', "clust2" ='#0099B4FF','clust3'='#42B540FF')))
#Z-score
tcga.pathway.exp_res1=t(apply(tcga.pathway.exp_res[,-c(1,2)],1,function(x){return(mosaic::zscore(x))}))
tcga.pathway.exp_res=cbind.data.frame(tcga.pathway.exp_res[,c(1,2)],
                                      tcga.pathway.exp_res1)
#第一个通路
pathway.gene1=cellage.pathway[cellage.pathway$ont %in% sig_pathway[1],]$gene
mat1=tcga.pathway.exp_res[tcga.pathway.exp_res$gene %in% pathway.gene1,]
rownames(mat1)=mat1$lab
#列注释
sig_pathway[1]
pathway= data.frame(pathway=sig_pathway[1],gene=rownames(mat1))
col_anno <-HeatmapAnnotation(pathway=pathway$pathway)
names(col_anno)=sig_pathway[1]
p1= Heatmap(matrix = as.matrix(t(mat1[,rownames(sam_anno)])),
            name="Row-Zscore Expression",
            cluster_columns = T,
            cluster_rows = F,
            column_title =sig_pathway[1],
            left_annotation= sam_ann,
            col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
            show_column_names = T,
            show_row_names = F,
            show_heatmap_legend = T,
            top_annotation = col_anno, 
            show_row_dend = T)
p1

for (i in 2:length(sig_pathway)){
  print(i)
  pathway.gene1=cellage.pathway[cellage.pathway$ont %in% sig_pathway[1],]$gene
  mat1=tcga.pathway.exp_res[tcga.pathway.exp_res$gene %in% pathway.gene1,]
  rownames(mat1)=mat1$lab
  #列注释
  #sig_pathway[i]
  pathway= data.frame(pathway=sig_pathway[i],gene=rownames(mat1))
  col_anno <-HeatmapAnnotation(pathway=pathway$pathway)
  names(col_anno)=sig_pathway[i]
  p1=p1+ Heatmap(matrix = as.matrix(t(mat1[,rownames(sam_anno)])),
              name="Row-Zscore Expression",
              cluster_columns = T,
              cluster_rows = F,
              column_title =sig_pathway[i],
              #left_annotation= sam_ann,
              col = colorRamp2(c(-3,0,3), c("blue", "white", "red")),
              show_column_names = T,
              show_row_names = F,
              show_heatmap_legend = T,
              top_annotation = col_anno, 
              show_row_dend = T)
}

p1
pdf('results/ana/05.pathway.exp/Fig5.pdf',height = 7,width = 30)
print(p1)
dev.off()
sig_boxplot<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab('')+labs(color=leg)
  return(pp)
}
head(tcga.cellage.score)
tcga.cellage.score.clust<-data.frame(tcga.cellage.score[tcga.subtype$Samples,],
                                     Cluster=tcga.subtype$Cluster)
#第一个通路
sig_pathway
fig5b1<-sig_boxplot(dat = tcga.cellage.score.clust[,c('Cluster',sig_pathway[1])],
                    leg = 'Cluster',
                    ylab = paste0(sig_pathway[1],' Score'),
                    palette = ggsci::pal_lancet()(9)[c(1,4,3)])
fig5b1
fig5b2<-sig_boxplot(dat = tcga.cellage.score.clust[,c('Cluster',sig_pathway[2])],
                    leg = 'Cluster',
                    ylab = paste0(sig_pathway[2],' Score'),
                    palette = ggsci::pal_lancet()(9)[c(1,4,3)])
fig5b2
fig5b3<-sig_boxplot(dat = tcga.cellage.score.clust[,c('Cluster',sig_pathway[3])],
                    leg = 'Cluster',
                    ylab = paste0(sig_pathway[3],' Score'),
                    palette = ggsci::pal_lancet()(9)[c(1,4,3)])
fig5b3

#7、临床表型的差异分析####
dir.create('results/ana/06.subtype.cli')
tcga_subtype.cli=merge(tcga.subtype,tcga_cli,by='Samples')
rownames(tcga_subtype.cli)=tcga_subtype.cli$Samples
tcga_subtype.cli$Event=ifelse(tcga_subtype.cli$OS==0,'Alive','Dead')
fivenum(na.omit(tcga_subtype.cli$Age))
tcga_subtype.cli$Age1=ifelse(tcga_subtype.cli$Age>55,'>55',ifelse(tcga_subtype.cli$Age<=55,'<=55',NA))
#
table(tcga_subtype.cli$T.Stage)
tcga_subtype.cli$T.Stage=gsub('[abcd]','',tcga_subtype.cli$T.Stage)
# 
# tcga_subtype.cli$T.Stage[tcga_subtype.cli$T.Stage=='T1'|tcga_subtype.cli$T.Stage=='T2']<-'T1+T2'
# tcga_subtype.cli$T.Stage[tcga_subtype.cli$T.Stage=='T3'|tcga_subtype.cli$T.Stage=='T4']<-'T3+T4'

table(tcga_subtype.cli$N.Stage)
tcga_subtype.cli$N.Stage=gsub('[abcd]','',tcga_subtype.cli$N.Stage)
tcga_subtype.cli$N.Stage=gsub('mi','',tcga_subtype.cli$N.Stage)
tcga_subtype.cli$N.Stage=gsub(' \\(i-\\)','',tcga_subtype.cli$N.Stage)
tcga_subtype.cli$N.Stage=gsub(' \\(i\\+\\)','',tcga_subtype.cli$N.Stage)

#tcga_subtype.cli$N.Stage[tcga_subtype.cli$N.Stage=='N1'|tcga_subtype.cli$N.Stage=='N2'|tcga_subtype.cli$N.Stage=='N3']<-'N1+M2+N3'

table(tcga_subtype.cli$M.Stage)
tcga_subtype.cli$M.Stage[tcga_subtype.cli$M.Stage=='MX'|tcga_subtype.cli$M.Stage=='']<-NA

table(tcga_subtype.cli$Stage)
tcga_subtype.cli$Stage=gsub('Stage ','',tcga_subtype.cli$Stage)
tcga_subtype.cli$Stage=gsub('[ABC]','',tcga_subtype.cli$Stage)
tcga_subtype.cli$Stage[tcga_subtype.cli$Stage=='']<-NA

# tcga_subtype.cli$Stage[tcga_subtype.cli$Stage=='I'|tcga_subtype.cli$Stage=='II']<-'I+II'
# tcga_subtype.cli$Stage[tcga_subtype.cli$Stage=='III'|tcga_subtype.cli$Stage=='IV']<-'III+IV'


table(tcga_subtype.cli$Gender)

table(tcga_subtype.cli$Cluster)

write.table(tcga_subtype.cli,'results/ana/06.subtype.cli/tcga_subtype.cli.txt',quote = F,row.names = F,sep='\t')


pie_compare_plot<-function(dat,gname,group_cols){
  library(dplyr)
  g_n<-as.character(unique(dat[,gname]))
  vname <- setdiff(colnames(dat), gname)
  pie.gname=data.frame()
  fisher.p <- c()
  for (i in vname) {
    tmp <- table(dat[,gname], dat[,i])
    p <- round(chisq.test(tmp)$p.value,4)
    names(p) <- i
    fisher.p <- c(fisher.p, p)
    pie.dat <- 
      tmp %>% as.data.frame() %>% group_by(Var1) %>% mutate(Pct = Freq/sum(Freq)) %>% as.data.frame()
    pie.gname=rbind.data.frame(pie.gname,pie.dat)
  }
  #plot
  vname_col<-list()
  for(i in 1:length(vname)){
    col_i<-length(as.character(unique(na.omit(dat[,vname[i]]))))
    vname_col[[i]] <- ggplot2::alpha(group_cols[i], (1:col_i)/col_i)
  }
  
  
  #names(vname_col)=vname
  col_num<-ncol(dat)
  row_num<-length(g_n)
  row_nums=1+2*row_num+2
  #c(1:col_num,2*row_num*col_num+1)
  #第一行
  nums1=c()
  for (i in 1:col_num){
    nums1=c(nums1,rep(i,3))
  }
  #两行
  nums2=c()
  for (j in 1:row_num){
    nums21=c()
    for (i in (col_num*j+1):((1+j)*col_num)){
      nums21=c(nums21,rep(i,3))
    }
    nums2=c(nums2,nums21,nums21)
  }
  
  #倒数第二行
  nums3=c()
  #(1+row_num)*col_num+1,(2+row_num)*col_num
  for (i in (((1+row_num)*col_num+1):((2+row_num)*col_num))){
    nums3=c(nums3,rep(i,3))
  }
  #最后一行
  nums4=c(rep(((2+row_num)*col_num)+1,col_num*3))
  nums=c(nums1,nums2,nums3,nums4)
  showLayout <- F 
  
  layout(matrix(nums,byrow = T,nrow = row_nums))
  if(showLayout) {
    layout.show(n = ((2+row_num)*col_num)+1) # 直观展示画布分布
  }
  par(bty="n", mgp = c(0,0,0), mar = c(0,0,0,0), lwd = 2) 
  plot(1,1,
       xlab = "",xaxt = "n", # 不显示x坐标轴
       ylab = "",yaxt = "n") # 不显示y坐标轴
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") # 背景涂黑
  text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
       (par("usr")[3]+par("usr")[4])/2,
       gname,cex = 2, col = "black") # 显示图标题
  #标题
  for (i in 1:length(vname)){
    plot(1,1,
         xlab = "",xaxt = "n", 
         ylab = "",yaxt = "n") 
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
    text((par("usr")[1]+par("usr")[2])/2, 
         (par("usr")[3]+par("usr")[4])/2,
         vname[i],cex = 2, col = "black") 
  }
  #抬头和扇形图
  for (i in 1:length(g_n)){
    plot(1,1,
         xlab = "",xaxt = "n", # 不显示x坐标轴
         ylab = "",yaxt = "n") # 不显示y坐标轴
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") # 背景涂黑
    text((par("usr")[1]+par("usr")[2])/2,
         (par("usr")[3]+par("usr")[4])/2,
         paste0(g_n[i],"\n(n = ",as.numeric(table(dat[,gname])[i]),")"),cex = 2, col = "black") 
    for (j in 1:length(vname)){
      aa=as.character(unique(dat[,vname[j]]))
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i],]
      pie.gname1=pie.gname[pie.gname$Var1==g_n[i] & pie.gname$Var2 %in% aa,]
      pie(pie.gname1$Pct, 
          col = vname_col[[j]], 
          border = "white",  
          radius = 1, 
          labels = NA,
          init.angle = 90)
      symbols(0,0, circles = .55, inches = FALSE, col = "white", bg = "white", lty = 0, add = TRUE)
    }
    
  }
  
  plot(1,1,
       xlab = "",xaxt = "n",
       ylab = "",yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white") 
  text((par("usr")[1]+par("usr")[2])/2,
       (par("usr")[3]+par("usr")[4])/2,
       'chisq.test',cex = 2, col = "black") 
  for(i in vname){
    plot(1,1,col = "white",
         xlab = "",xaxt = "n", # 不显示x坐标轴
         ylab = "",yaxt = "n") # 不显示y坐标轴
    text((par("usr")[1]+par("usr")[2])/2, # 用par("usr")来获取该画布的绝对位置
         (par("usr")[3]+par("usr")[4])/2,
         paste0("p = ",fisher.p[i]),cex = 1.5, col = "black") # 显示图标题
    abline(h = par("usr")[3], col = "black")
  }
  plot(0,0,col = "white",
       xlab = "",xaxt = "n", # 不显示x坐标轴
       ylab = "",yaxt = "n") # 不显示y坐标轴
  leg_nam=c()
  for (ii in 1:length(vname)){
    aa=as.character(unique(dat[,vname[ii]]))
    pie.gname1=pie.gname[pie.gname$Var2 %in% aa,]
    leg_nam=c(leg_nam,as.character(unique(pie.gname1$Var2)))
  }
  vname_col_name=unlist(vname_col)
  
  legend("topleft",ncol=ceiling(length(leg_nam)/2)+1,
         legend = leg_nam,
         fill = vname_col_name,
         border = NA,
         bty = "n", 
         cex = 1.2,
         x.intersp = 0.05,
         y.intersp = 1,
         text.width = 0.1, 
         horiz = F)
}


pdf('results/ana/06.subtype.cli/Fig6.pdf',height = 7,width = 15)
pie_compare_plot(dat = tcga_subtype.cli[,c("Cluster","T.Stage","N.Stage","Stage","Age1",'Event')],
                 gname = 'Cluster',
                 group_cols = c(ggsci::pal_lancet()(9)))
dev.off()

#8、突变的差异####
dir.create('results/ana/07.subtype.cnv/clust1',recursive = T)
dir.create('results/ana/07.subtype.cnv/clust2',recursive = T)
dir.create('results/ana/07.subtype.cnv/clust3',recursive = T)

metadata <- read.delim('/pub1/data/mg_projects/TCGA/Matrix/cnvs/all.sample.txt',sep='\t',header = F,check.names = F)
cnv_dir='/pub1/data/mg_projects/TCGA/Matrix/cnvs/'
tur_name='TCGA-BRCA'
metadata=metadata[grep(pattern = tur_name,metadata$V2),c(2,3,4)]
head(metadata)

metadata=data.frame(GDC_Aliquot=stringr::str_split_fixed(metadata$V2,'/',6)[,6],
                    Samples=paste0(metadata$V3,'-',ifelse(metadata$V4==1,'01',metadata$V4)))
head(metadata)
metadata=metadata[grep('nocnv_grch38.seg.v2.txt$',metadata$GDC_Aliquot),]
cnv_df=data.frame()
for (i in 1:nrow(metadata)){
  cnv_df1=data.table::fread(file=paste0(cnv_dir,'/',tur_name,'/',metadata$GDC_Aliquot[i]),
                            data.table = F)
  cnv_df1$Samples=metadata$Samples[i]
  cnv_df=rbind(cnv_df,cnv_df1)
  
}

head(cnv_df)
cnv_df <- cnv_df[,c("Samples", "Chromosome","Start", "End", "Num_Probes", "Segment_Mean")]
head(cnv_df)

table(tcga.subtype$Cluster)
colnames(cnv_df)[1]='Sample'
c1_sam<-tcga.subtype[tcga.subtype$Cluster=='clust1',"Samples"]
c2_sam<-tcga.subtype[tcga.subtype$Cluster=='clust2',"Samples"]
c3_sam<-tcga.subtype[tcga.subtype$Cluster=='clust3',"Samples"]

cnv_df=unique(cnv_df)
cnv_df1=cnv_df[cnv_df$Sample %in% c1_sam,]
cnv_df2=cnv_df[cnv_df$Sample %in% c2_sam,]
cnv_df3=cnv_df[cnv_df$Sample %in% c3_sam,]

cnv_df1=na.omit(cnv_df1)
cnv_df2=na.omit(cnv_df2)
cnv_df3=na.omit(cnv_df3)

#过滤
write.table(cnv_df1, file = "results/ana/07.subtype.cnv/C1_segment_file.txt", sep = "\t", row.names = F, quote = F)
write.table(cnv_df2, file = "results/ana/07.subtype.cnv/C2_segment_file.txt", sep = "\t", row.names = F, quote = F)
write.table(cnv_df3, file = "results/ana/07.subtype.cnv/C3_segment_file.txt", sep = "\t", row.names = F, quote = F)




#绘图
library(TCGAbiolinks)
library(SummarizedExperiment)
chrom_extract <- function(BSgenome.hg  = NULL) {
  if (is.null(BSgenome.hg )) stop("NULL object !", call. = FALSE)
  obj <- list(species = GenomeInfoDb::organism(BSgenome.hg), genomebuild = BSgenome::providerVersion(BSgenome.hg))
  df <- data.frame(chrom = BSgenome::seqnames(BSgenome.hg), chrN = seq_along(BSgenome::seqnames(BSgenome.hg)), chr.length = GenomeInfoDb::seqlengths(BSgenome.hg), stringsAsFactors = FALSE)
  df <- df[1:24,]
  df$chr.length.sum <- cumsum(as.numeric(df$chr.length))
  df$chr.length.cumsum <- c(0, df$chr.length.sum[-nrow(df)])
  df$middle.chr <- round(diff(c(0, df$chr.length.sum)) /2)
  df$middle.chr.genome <- df$middle.chr + df$chr.length.cumsum
  obj$chromosomes <- df
  obj$chrom2chr <- sapply(obj$chromosomes$chrom, function(k) { obj$chromosomes$chrN[obj$chromosomes$chrom == k]}, simplify = FALSE)
  obj$chr2chrom <- sapply(obj$chromosomes$chrN, function(k) { obj$chromosomes$chrom[obj$chromosomes$chrN == k]}, simplify = FALSE)
  names(obj$chr2chrom) <- obj$chromosomes$chrN
  obj$genome.length <- sum(as.numeric(obj$chromosomes$chr.length), na.rm = TRUE)
  return(obj)
}

# Extract a chromosomes reference loci
BSgenome.hg = "BSgenome.Hsapiens.UCSC.hg38"
BSg.obj <- getExportedValue(BSgenome.hg, BSgenome.hg)
genome.version <- BSgenome::providerVersion(BSg.obj)
chrom <- chrom_extract(BSg.obj)
#clsut1
C1_scores <- read.table("results/ana/07.subtype.cnv/clust1/scores.gistic", sep="\t",header=T,stringsAsFactors = F)
C1_scores=C1_scores[C1_scores$Type %in% c('Amp','Del'),]
#clust2
C2_scores <- read.table("results/ana/07.subtype.cnv/clust2/scores.gistic", sep="\t",header=T,stringsAsFactors = F)
C2_scores=C2_scores[C2_scores$Type %in% c('Amp','Del'),]
#clust3
C3_scores <- read.table("results/ana/07.subtype.cnv/clust3/scores.gistic", sep="\t",header=T,stringsAsFactors = F)
C3_scores=C3_scores[C3_scores$Type %in% c('Amp','Del'),]


cnv_freq_plot<-function(dat,yname,ylab,col,chrom,ylim){
  dat[dat$Chromosome==23,"Chromosome"]="X"
  dat[dat$Chromosome==24,"Chromosome"]="Y"
  chrID <- unname(unlist(chrom$chrom2chr[as.character(paste0("chr",dat$Chromosome))]))
  dat$Start.geno <- dat$Start + chrom$chromosomes$chr.length.cumsum[chrID]
  dat$End.geno <- dat$End + chrom$chromosomes$chr.length.cumsum[chrID]
  dat.amp <- dat[dat$Type=="Amp",]
  dat.del <- dat[dat$Type=="Del",]
  ink <- chrom$chromosomes$chrN %in% chrID
  
  if(yname=='G.score'){
    m.pos <- c(ylim[1]+1,ylim[2]-1)
    m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
  }else{
    m.pos <- c(ylim[1]+10,ylim[2]-10)
    m.mod <- -(chrom$chromosomes$chrN[ink] %% 2) +2
  }
  plot(dat.amp$Start.geno, dat.amp[,yname],
       pch = ".", type='h',cex = 2, xaxs = "i", yaxs = "i",
       xlim = c(0,chrom$genome.length), 
       ylim = ylim,
       main = "", cex.main = 2, 
       ylab = ylab, xlab = NA,
       cex.lab = 1, col = adjustcolor(col[1], alpha.f = .8),
       xaxt = "n", lwd = 2, las=1) 
  lines(dat.del$Start.geno, dat.del[,yname], type='h', lwd = 2, col = adjustcolor(col[2], alpha.f = .8))
  yrange = abs(diff(ylim))+2
  try(text(x = chrom$chromosomes$middle.chr.geno[ink], 
           y = m.pos[m.mod], 
           labels = stringr::str_sub(chrom$chromosomes$chrom[ink],
                                     start = 4), cex = 1))
  abline(h = 0.0, col = 1, lwd = 1, lty = 3)
  abline(v = c(0,chrom$chromosomes$chr.length.sum), col = 1, lty = 3, lwd = 1)
  
  col1 <- adjustcolor(col[1], alpha.f = .8)
  col2 <- adjustcolor(col[2], alpha.f = .8)
  legend("topleft", c("gain","loss"), cex=0.6, bty="n", fill=c(col1,col2))
}

#freq
C1_scores$frequency=ifelse(C1_scores$Type=='Amp',C1_scores$frequency *100,C1_scores$frequency * (-100))
C2_scores$frequency=ifelse(C2_scores$Type=='Amp',C2_scores$frequency *100,C2_scores$frequency * (-100))
C3_scores$frequency=ifelse(C3_scores$Type=='Amp',C3_scores$frequency *100,C3_scores$frequency * (-100))

pdf("results/ana/07.subtype.cnv/C1_copy_number_frequency.pdf",width = 10,height = 5)
cnv_freq_plot(dat = C1_scores,yname = 'frequency',
              ylab='Percentage/Frequency',
              ylim = c(-100,100),chrom = chrom,
              col=ggsci::pal_npg("nrc")(10))
dev.off()
pdf("results/ana/07.subtype.cnv/C2_copy_number_frequency.pdf",width = 10,height = 5)
cnv_freq_plot(dat = C2_scores,yname = 'frequency',
              ylab='Percentage/Frequency',
              ylim = c(-100,100),chrom = chrom,
              col=ggsci::pal_npg("nrc")(10))
dev.off()
pdf("results/ana/07.subtype.cnv/C3_copy_number_frequency.pdf",width = 10,height = 5)
cnv_freq_plot(dat = C3_scores,yname = 'frequency',
              ylab='Percentage/Frequency',
              ylim = c(-100,100),chrom = chrom,
              col=ggsci::pal_npg("nrc")(10))
dev.off()
#snv
tcga_maf=getTCGAMAFByCode('BRCA')

#查看
tcga.subtype1=tcga.subtype[,c("Samples","Cluster")]
colnames(tcga.subtype1)[1]='Tumor_Sample_Barcode'
tcga.subtype1$Tumor_Sample_Barcode=substr(tcga.subtype1$Tumor_Sample_Barcode,1,12)

tcga.maf1<-read.maf('origin_datas/TCGA/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf',
                    isTCGA=T)
tcga.maf1=subsetMaf(tcga.maf1,tsb=intersect(tcga.maf1@data$Tumor_Sample_Barcode,tcga.subtype1$Tumor_Sample_Barcode))

write.table(tcga.subtype1[,c("Tumor_Sample_Barcode","Cluster")],'results/ana/07.subtype.cnv/tcga.subtype.txt',quote = F,row.names = F,sep='\t')

tcga.maf1<-read.maf(tcga.maf1@data,
                    isTCGA=T,clinicalData = 'results/ana/07.subtype.cnv/tcga.subtype.txt')

tcga.maf1@clinical.data
ggsci::pal_lancet()(9)[c(1,4,3)]
pdf('results/ana/07.subtype.cnv/Fig7B.pdf',width = 10,height = 5)
oncoplot(maf=tcga.maf1,clinicalFeatures = 'Cluster',
         top = 10,sortByAnnotation = T,
         annotationColor = list(Cluster=c(clust1='#00468BFF',clust2='#0099B4FF',clust3='#42B540FF')))
#colors =ggsci::pal_lancet()(9)[3:5])
dev.off()

#9、生物学特征####
dir.create('results/ana/08.bio.feat')

tcga_tmp_log2_T=tcga_tmp_log2[,tcga_cli$Samples]
dim(tcga_tmp_log2_T)
#CCP
library(ggpubr)
immu_ccp<-function(exp_data){
  #cell cycle progression(CCP,31个基因 ,PMID: 21310658)
  ccp.gene=c('FOXM1','ASPM','TK1','PRC1','CDC20','BUB1B','PBK','DTL','CDKN3','RRM2','ASF1B','CEP55','CDC2','DLGAP5','C18orf24','RAD51','KIF11','BIRC5','RAD54L','CENPM','KIAA0101','KIF20A','PTTG1','CDCA8','NUSAP1','PLK1','CDCA3','ORC6L','CENPF','TOP2A','MCM10')
  gs.score=ssGSEAScore_by_genes(exp_data,ccp.gene)
  rownames(gs.score)='CCP'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga_tmp_log2_T=tcga_tmp_log2[,tcga.subtype$Samples]
tcga.ccp=immu_ccp(tcga_tmp_log2_T)
tcga.ccp=crbind2DataFrame(tcga.ccp)
tcga.ccp.subtype=merge(data.frame(Samples=rownames(tcga.ccp),CCP=tcga.ccp$CCP),
                       tcga.subtype,by='Samples')
tcga.ccp.subtype=tcga.ccp.subtype[order(tcga.ccp.subtype$Cluster),]

fig8a<-sig_boxplot(dat = tcga.ccp.subtype[,c("Cluster","CCP")],
                   leg = 'Cluster',ylab = 'CCP Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8a
write.table(tcga.ccp.subtype,'results/ana/08.bio.feat/tcga.ccp.subtype.txt',quote = F,sep='\t',row.names = F)
#Cell cycle G1/S(https://www.kegg.jp/pathway/hsa04110)
immu_g1.s.cell.cycle<-function(exp_data){
  #G1/S cell cycle(27个)
  gene=c('MYC','MAX','ZBTB17','CDKN2B','CCND1','CCND2','CCND3','CDK4','CDK6','RB1','E2F1','E2F2','E2F3','CDKN2A','MDM2','TP53','CDKN1A','SKP1','CUL1','RBX1','SKP2','CDKN1B','CCNA2','CCNA1','CCNE1','CCNE2','CDK2')
  gs.score=ssGSEAScore_by_genes(exp_data,gene)
  rownames(gs.score)='G1/S Cell Cycle'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga.g1.s.cell.cycle=immu_g1.s.cell.cycle(tcga_tmp_log2_T)
tcga.g1.s.cell.cycle=crbind2DataFrame(tcga.g1.s.cell.cycle)
tcga.g1.s.cell.cycle.subtype=merge(data.frame(Samples=rownames(tcga.g1.s.cell.cycle),pathway=tcga.g1.s.cell.cycle$`G1/S Cell Cycle`),
                                   tcga.subtype,by='Samples')
tcga.g1.s.cell.cycle.subtype=tcga.g1.s.cell.cycle.subtype[order(tcga.g1.s.cell.cycle.subtype$Cluster),]

fig8b<-sig_boxplot(dat = tcga.g1.s.cell.cycle.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'G1/S Cell Cycle Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8b
write.table(tcga.g1.s.cell.cycle.subtype,'results/ana/08.bio.feat/tcga.g1.s.cell.cycle.subtype.txt',quote = F,sep='\t',row.names = F)
#

#G2M_CHECKPOINT
immu_g2m.checkpoint<-function(exp_data){
  #HALLMARK_G2M_CHECKPOINT h.all.v7.0.symbols.gmt
  gene=c('AURKA','CCNA2','TOP2A','CCNB2','CENPA','BIRC5','CDC20','PLK1','TTK','PRC1','NDC80','KIF11','NUSAP1','CKS2','KIF2C','MKI67','AURKB','TPX2','SMC4','BUB1','CENPF','RACGAP1','CENPE','AC027237.1','UBE2C','MCM6','MCM3','PTTG1','CDK1','KIF4A','ESPL1','MAD2L1','NEK2','KIF22','HMMR','KPNA2','CDKN3','CDC25A','H2AFX','CDC25B','PLK4','CDC6','CCNF','MCM5','LMNB1','E2F3','KIF15','CHEK1','UBE2S','NSD2','HMGB3','DBF4','TACC3','MCM2','CDKN2C','CDKN1B','FANCC','NASP','STAG1','GINS2','FBXO5','POLQ','EZH2','RAD21','STMN1','SUV39H1','PRIM2','E2F1','CHAF1A','NOLC1','GSPT1','BUB3','SMC1A','ILF3','CDC7','INCENP','CKS1B','EXO1','H2AFZ','TFDP1','CCND1','KPNB1','JPT1','LBR','HUS1','KIF20B','TOP1','PDS5B','SRSF1','STIL','ABL1','DTYMK','CDC27','BARD1','ATF5','CDC45','ODC1','XPO1','SFPQ','TMPO','PML','BRCA2','CTCF','KNL1','KMT5A','SLC38A1','TRA2B','MYBL2','TROAP','TENT4A','CUL3','MAPK14','HIST1H2BK','MYC','AMD1','CBX1','CHMP1A','DKC1','YTHDC1','CCNT1','TGFB1','ATRX','LIG3','NUP50','SLC7A5','RBL1','NUMA1','RAD54L','EFNA5','PRPF4B','UCK2','ARID4A','CUL1','UPF1','DR1','MNAT1','SMC2','RBM14','RPA2','SQLE','ORC6','CDK4','POLE','RASAL2','HOXC10','RPS6KA5','CUL4A','SLC7A1','FOXN3','HMGA1','AC091021.1','TRAIP','PRMT5','CUL5','DDX39A','MARCKS','PBK','ORC5','SAP30','KATNA1','HNRNPD','POLA2','HIRA','HIF1A','SYNCRIP','TLE3','NCL','RAD23B','E2F2','EGF','HMGN2','SRSF10','SNRPD1','CASP8AP2','SMARCC1','SLC12A2','NOTCH2','TNPO2','SMAD3','MAP3K20','HSPA8','G3BP1','PTTG3P','DMD','MEIS1','HNRNPU','SRSF2','MT2A','NUP98','EWSR1','KIF5B','MTF2','E2F4','BCL3','PURA','MEIS2','PAFAH1B1','WRN','H2AFV','ODF2')
  gs.score=ssGSEAScore_by_genes(exp_data,gene)
  rownames(gs.score)='G2M Checkpoint'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga.g2m.checkpoint=immu_g2m.checkpoint(tcga_tmp_log2_T)
tcga.g2m.checkpoint=crbind2DataFrame(tcga.g2m.checkpoint)
tcga.g2m.checkpoint.subtype=merge(data.frame(Samples=rownames(tcga.g2m.checkpoint),pathway=tcga.g2m.checkpoint$`G2M Checkpoint`),
                                  tcga.subtype,by='Samples')
tcga.g2m.checkpoint.subtype=tcga.g2m.checkpoint.subtype[order(tcga.g2m.checkpoint.subtype$Cluster),]

fig8c<-sig_boxplot(dat = tcga.g2m.checkpoint.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'G2M Checkpoint Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8c
write.table(tcga.g2m.checkpoint.subtype,'results/ana/08.bio.feat/tcga.g2m.checkpoint.subtype.txt',quote = F,sep='\t',row.names = F)
#端粒酶延长
#REACTOME_TELOMERE_EXTENSION_BY_TELOMERAS_MERASE
immu_telomere<-function(exp_data){
  #HALLMARK_G2M_CHECKPOINT h.all.v7.0.symbols.gmt
  gene=c('ACD','ANKRD28','CCNA1','CCNA2','CDK2','DKC1','GAR1','NHP2','NOP10','PIF1','POT1','PPP6C','PPP6R3','RTEL1','RUVBL1','RUVBL2','SHQ1','TERF1','TERF2','TERF2IP','TERT','TINF2','WRAP53')
  gs.score=ssGSEAScore_by_genes(exp_data,gene)
  rownames(gs.score)='Telomere Extension By Telomerase '
  gs.score=t(gs.score)
  return(gs.score)
}
tcga.telomere=immu_telomere(tcga_tmp_log2_T)
tcga.telomere=crbind2DataFrame(tcga.telomere)
tcga.telomere.subtype=merge(data.frame(Samples=rownames(tcga.telomere),pathway=tcga.telomere$`Telomere Extension By Telomerase `),
                            tcga.subtype,by='Samples')
tcga.telomere.subtype=tcga.telomere.subtype[order(tcga.telomere.subtype$Cluster),]

fig8d<-sig_boxplot(dat = tcga.telomere.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'Telomere Extension By Telomerase Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8d
write.table(tcga.telomere.subtype,'results/ana/08.bio.feat/tcga.telomere.subtype.txt',quote = F,sep='\t',row.names = F)

#EMT
immu_EMT<-function(exp_data){
  #HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION :h.all.v7.0.symbols.gmt
  EMT.gene=c('COL3A1','COL5A2','COL5A1','FBN1','COL1A1','FN1','COL6A3','SERPINE1','COL1A2','COL4A1','COL4A2','VCAN','IGFBP3','TGFBI','SPARC','LUM','LAMC1','LOX','LAMC2','CCN2','TAGLN','COL7A1','LOXL2','COL6A2','ITGAV','THBS2','COL16A1','NNMT','TPM1','CDH2','MMP2','COL11A1','THBS1','FAP','BGN','SERPINH1','FSTL1','POSTN','THY1','SPP1','TNC','TFPI2','NID2','ITGB5','MMP3','VIM','LOXL1','FBLN5','COL12A1','ELN','CDH11','COMP','SPOCK1','BMP1','IL32','LAMA3','TIMP1','QSOX1','TIMP3','VCAM1','CCN1','EDIL3','CALD1','MAGEE1','FBLN1','SGCB','ECM1','LAMA2','FSTL3','TPM2','INHBA','DAB2','EMP3','BASP1','ITGA5','MGP','VEGFA','CXCL1','WNT5A','SDC1','PLOD2','PCOLCE','GREM1','ITGB1','COL5A3','RHOB','HTRA1','FGF2','SNTB1','GADD45A','MEST','LRRC15','TNFRSF11B','CD59','ACTA2','EFEMP2','MATN2','PCOLCE2','SERPINE2','GPC1','ABI3BP','FUCA1','SLIT3','LAMA1','PMEPA1','COL8A2','FBN2','IGFBP2','PFN2','SDC4','CD44','GADD45B','CXCL8','GLIPR1','ANPEP','P3H1','VEGFC','MMP14','SGCD','PLOD1','MATN3','MYL9','SLC6A8','CALU','PRRX1','TNFRSF12A','FMOD','ID2','GEM','PLAUR','MYLK','TGFB1','SFRP1','PLOD3','IL6','APLP1','FBLN2','MSX1','PTX3','FZD8','JUN','FERMT2','DKK1','SNAI2','DST','TPM4','DCN','GJA1','PMP22','IGFBP4','COPA','LRP1','ITGA2','FLNA','MFAP5','PTHLH','TGFBR3','SFRP4','LGALS1','RGS4','CDH6','SAT1','NT5E','DPYSL3','PPIB','TGM2','SGCG','ITGB3','PDLIM4','CTHRC1','ECM2','CRLF1','AREG','IL15','MCM7','GAS1','PRSS2','CADM1','OXTR','SCG2','CXCL6','MMP1','TNFAIP3','CAPG','CAP2','MXRA5','FOXC2','NTM','ENO2','FAS','BDNF','ADAM12','PVR','CXCL12','PDGFRB','SLIT2','NOTCH2','COLGALT1','GPX7','WIPF1')
  gs.score=ssGSEAScore_by_genes(exp_data,EMT.gene)
  rownames(gs.score)='EMT'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga.EMT=immu_EMT(tcga_tmp_log2_T)
tcga.EMT=crbind2DataFrame(tcga.EMT)
tcga.EMT.subtype=merge(data.frame(Samples=rownames(tcga.EMT),EMT=tcga.EMT$EMT),
                       tcga.subtype,by='Samples')
write.table(tcga.EMT.subtype,'results/ana/08.bio.feat/tcga.EMT.subtype.txt',quote = F,row.names = F,sep='\t')
tcga.EMT.subtype=tcga.EMT.subtype[order(tcga.EMT.subtype$Cluster),]

fig8e<-sig_boxplot(dat = tcga.EMT.subtype[,c("Cluster","EMT")],
                   leg = 'Cluster',ylab = 'EMT Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8e
#缺氧
immu_hypoxia<-function(exp_data){
  #缺氧型基因 HALLMARK_HYPOXIA
  hypoxia.gene=c('PGK1','PDK1','GBE1','PFKL','ALDOA','ENO2','PGM1','NDRG1','HK2','ALDOC','GPI','MXI1','SLC2A1','P4HA1','ADM','P4HA2','ENO1','PFKP','AK4','FAM162A','PFKFB3','VEGFA','BNIP3L','TPI1','ERO1A','KDM3A','CCNG2','LDHA','GYS1','GAPDH','BHLHE40','ANGPTL4','JUN','SERPINE1','LOX','GCK','PPFIA4','MAFF','DDIT4','SLC2A3','IGFBP3','NFIL3','FOS','RBPJ','HK1','CITED2','ISG20','GALK1','WSB1','PYGM','STC1','ZNF292','BTG1','PLIN2','CSRP2','VLDLR','JMJD6','EXT1','F3','PDK3','ANKZF1','UGP2','ALDOB','STC2','ERRFI1','ENO3','PNRC1','HMOX1','PGF','GAPDHS','CHST2','TMEM45A','BCAN','ATF3','CAV1','AMPD3','GPC3','NDST1','IRS2','SAP30','GAA','SDC4','STBD1','IER3','PKLR','IGFBP1','PLAUR','CAVIN3','CCN5','LARGE1','NOCT','S100A4','RRAGD','ZFP36','EGFR','EDN2','IDS','CDKN1A','RORA','DUSP1','MIF','PPP1R3C','DPYSL4','KDELR3','DTNA','ADORA2B','HS3ST1','CAVIN1','NR3C1','KLF6','GPC4','CCN1','TNFAIP3','CA12','HEXA','BGN','PPP1R15A','PGM2','PIM1','PRDX5','NAGK','CDKN1B','BRS3','TKTL1','MT1E','ATP7A','MT2A','SDC3','TIPARP','PKP1','ANXA2','PGAM2','DDIT3','PRKCA','SLC37A4','CXCR4','EFNA3','CP','KLF7','CCN2','CHST3','TPD52','LXN','B4GALNT2','PPARGC1A','BCL2','GCNT2','HAS1','KLHL24','SCARB1','SLC25A1','SDC2','CASP6','VHL','FOXO3','PDGFB','B3GALT6','SLC2A5','SRPX','EFNA1','GLRX','ACKR3','PAM','TGFBI','DCN','SIAH2','PLAC8','FBP1','TPST2','PHKG1','MYH9','CDKN1C','GRHPR','PCK1','INHA','HSPA5','NDST2','NEDD4L','TPBG','XPNPEP1','IL6','SLC6A6','MAP3K1','LDHC','AKAP12','TES','KIF5A','LALBA','COL5A1','GPC1','HDLBP','ILVBL','NCAN','TGM2','ETS1','HOXB9','SELENBP1','FOSL2','SULT2B1','TGFB3')
  gs.score=ssGSEAScore_by_genes(exp_data,hypoxia.gene)
  rownames(gs.score)='hypoxia.score'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga_hypoxia<-immu_hypoxia(exp_data = tcga_tmp_log2_T)
tcga_hypoxia=crbind2DataFrame(tcga_hypoxia)

tcga_hypoxia.subtype=merge(data.frame(Samples=rownames(tcga_hypoxia),pathway=tcga_hypoxia$hypoxia.score),
                           tcga.subtype,by='Samples')
write.table(tcga_hypoxia.subtype,'results/ana/08.bio.feat/tcga_hypoxia.subtype.txt',quote = F,row.names = F,sep='\t')

tcga_hypoxia.subtype=tcga_hypoxia.subtype[order(tcga_hypoxia.subtype$Cluster),]

fig8f<-sig_boxplot(dat = tcga_hypoxia.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'Hypoxia Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])
fig8f

#血管生成
immu_Angiogenesis<-function(exp_data){
  #血管生成相关的基因集，文献来源：PMC3743050
  angiogenesis.sym=c('TEK','CDH5','DLL4','EGFL7','KDR','ROBO4','VWF','TMSB4XP2','EDIL3','PECAM1','CD34','ENG','ESM1','ADGRF5','SELE','PROM1','VEGFA','RAMP2','RGS5','AGGF1','ITGB3','ITGAV','ITGB5','ANGPTL4')
  gs.score=ssGSEAScore_by_genes(exp_data,angiogenesis.sym)
  rownames(gs.score)='angiogenesis.score'
  gs.score=t(gs.score)
  return(gs.score)
}
tcga_Angiogenesis<-immu_Angiogenesis(exp_data = tcga_tmp_log2_T)
tcga_Angiogenesis=crbind2DataFrame(tcga_Angiogenesis)

tcga_Angiogenesis.subtype=merge(data.frame(Samples=rownames(tcga_Angiogenesis),pathway=tcga_Angiogenesis$angiogenesis.score),
                                tcga.subtype,by='Samples')
write.table(tcga_Angiogenesis.subtype,'results/ana/08.bio.feat/tcga_Angiogenesis.subtype.txt',quote = F,row.names = F,sep='\t')

tcga_Angiogenesis.subtype=tcga_Angiogenesis.subtype[order(tcga_Angiogenesis.subtype$Cluster),]

fig8g<-sig_boxplot(dat = tcga_Angiogenesis.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'Angiogenesis Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])
fig8g


#10个肿瘤相关的通路
#tcga.onco_pathway=immu_oncogenic_pathways(tcga_tmp_log2_T)
#save(tcga.onco_pathway,file='tcga.onco_pathway.RData')
load('tcga.onco_pathway.RData')
head(gse.onco_pathway)
fig8h=mg_PlotMutiBoxplot(tcga.onco_pathway[tcga.subtype$Samples,],
                         tcga.subtype$Cluster
                         , test_method = 'kruskal.test'
                         , add = 'boxplot'
                         , legend.pos = 'top'
                         , group_cols =ggsci::pal_lancet()(9)[c(1,4,3)]
                         , ylab = 'Enrichment Score',xangle = 0)
fig8h
write.table(tcga.onco_pathway,'results/ana/08.bio.feat/tcga.onco_pathway.txt',quote = F,sep='\t',row.names = T)


#炎性因子
immu_inflamatory.response<-function(exp_data){
  #HALLMARK_INFLAMMATORY_RESPONSE :h.all.v7.0.symbols.gmt
  gene=c('CXCL10','CCL2','CCL5','FPR1','CCL20','IL1A','CXCL8','CCL7','CCL22','CXCL11','CCR7','EDN1','CD40','CXCL9','IL6','IL1B','TLR2','IL1R1','CD69','ICAM1','CCRL2','AQP9','EREG','C3AR1','GNA15','CMKLR1','PTGER4','LIF','IL15','NAMPT','OPRK1','ITGB8','PTAFR','ADM','PLAUR','NFKB1','INHBA','OSM','TNFSF10','TNFSF15','IFNGR2','ADGRE1','IL12B','CSF1','CXCL6','TNFRSF9','LYN','ACVR2A','LDLR','BDKRB1','HRH1','F3','BST2','PTGIR','CD55','CALCRL','CSF3','GPR132','IL4R','NLRP3','IL15RA','ADORA2B','GCH1','OLR1','PTGER2','CSF3R','MYC','RELA','TNFAIP6','IL7R','IL18','GABBR1','CD82','TNFSF9','NMUR1','IL2RB','TLR1','LPAR1','IRAK2','RIPK2','MMP14','P2RX7','SLC11A2','SELL','P2RY2','ABCA1','FFAR2','PROK2','GNAI3','TACR1','SLC7A1','CDKN1A','CYBB','TIMP1','HBEGF','SCARF1','EBI3','NFKBIA','SRI','SLC7A2','CCL17','TLR3','APLNR','OSMR','IL10RA','PSEN1','GPR183','ATP2B1','TNFRSF1B','BEST1','GPC3','SCN1B','ACVR1B','HPN','SEMA4D','KLF6','CD48','CXCR6','SLC1A2','GP1BA','TAPBP','RGS16','SLAMF1','LCK','HIF1A','AHR','NMI','RHOG','TPBG','NPFFR2','IFNAR1','ICOSLG','RASGRP1','IFITM1','KCNJ2','LY6E','IL18R1','IL10','KCNA3','HAS2','DCBLD2','LAMP3','VIP','CD70','RGS1','SLC31A1','ADRM1','KCNMB2','SERPINE1','MXD1','AXL','MEFV','PVR','CCL24','PDE4B','LCP2','PDPN','IRF7','MET','ATP2A2','SLC31A2','FZD5','ITGA5','SGMS2','MARCO','CD14','EIF2AK2','ROS1','ATP2C1','NDP','BTG2','MSR1','PTPRE','RNF144B','PCDH7','SPHK1','IL18RAP','RTP4','RAF1','CHST2','ITGB3','KIF1B','SELE','NOD2','C5AR1','EMP3','CLEC5A','TACR3','SLC4A4','MEP1A','SELENOS','LTA','PIK3R5','STAB1','IRF1','ICAM4','P2RX4','ABI1','CX3CL1','SLC28A2')
  gs.score=ssGSEAScore_by_genes(exp_data,gene)
  rownames(gs.score)='inflamatory.response'
  gs.score=t(gs.score)
  return(gs.score)
}

tcga.inflamatory.response=immu_inflamatory.response(tcga_tmp_log2_T)
tcga.inflamatory.response=crbind2DataFrame(tcga.inflamatory.response)
tcga.inflamatory.response.subtype=merge(data.frame(Samples=rownames(tcga.inflamatory.response),pathway=tcga.inflamatory.response$inflamatory.response),
                                        tcga.subtype,by='Samples')
tcga.inflamatory.response.subtype=tcga.inflamatory.response.subtype[order(tcga.inflamatory.response.subtype$Cluster),]

fig8I<-sig_boxplot(dat = tcga.inflamatory.response.subtype[,c("Cluster","pathway")],
                   leg = 'Cluster',ylab = 'Inflamatory Response Score',
                   palette = ggsci::pal_lancet()(9)[c(1,4,3)])

fig8I
write.table(tcga.inflamatory.response.subtype,'results/ana/08.bio.feat/tcga.inflamatory.response.subtype.txt',quote = F,sep='\t',row.names = F)

fig8af<-mg_merge_plot(fig8a,fig8b,fig8c,fig8d,fig8e,fig8f,nrow = 1,ncol = 6,labels = LETTERS[1:6],common.legend = T)
fig8gi<-mg_merge_plot(fig8g,fig8I,fig8h,nrow = 1,ncol = 3,widths = c(1,1,3),labels = LETTERS[7:9],common.legend = T)
fig8<-mg_merge_plot(fig8af,fig8gi,nrow = 2,ncol = 1,common.legend = T)
ggsave('results/ana/08.bio.feat/Fig8.pdf',fig8,height = 9,width = 15)
#10、免疫特征#####
dir.create('results/ana/09.subtype.immnu')
#ESTIMATE
#tcga_esti<-immu_estimate(exp = tcga_tmp_log2_T,platform='illumina',isTCGA=T)
#save(tcga_esti,file='tcga_esti.RData')
load('tcga_esti.RData')
fig9a<-mg_PlotMutiBoxplot(data = tcga_esti[tcga.subtype$Samples,1:3],
                          group = tcga.subtype$Cluster,
                          group_cols = ggsci::pal_lancet()(9)[c(1,4,3)],
                          test_method = 'kruskal.test',
                          xangle = 30,
                          add = 'boxplot', 
                          legend.pos = 'top',ylab = 'Immnune score')
fig9a
write.table(tcga_esti,'results/ana/09.subtype.immnu/tcga_esti.txt',quote = F,sep='\t',row.names = T)

#tcga.cier<-immu_CIBERSORT(exp_data = tcga_tmp_log2_T)
#save(tcga.cier,file='tcga.cier.RData')
# load('tcga.cier.RData')
# tcga.cier[1:4,1:4]
# fig9b<-mg_PlotMutiBoxplot(tcga.cier[tcga.subtype$Samples,1:22],
#                           tcga.subtype$Cluster
#                           , test_method = 'kruskal.test'
#                           , add = 'boxplot'
#                           , legend.pos = 'top'
#                           , group_cols = ggsci::pal_lancet()(9)[c(1,4,3)]
#                           , ylab = 'Enrichment Score',xangle = 30)
# fig9b
#write.table(tcga.cier[,1:22],'results/ana/08.bio.feat/tcga.cier.txt',quote = F,sep='\t',row.names = T)
tcga.mcp<-immu_MCPcounter(exp = tcga_tmp_log2_T,isTCGA = T)
head(tcga.mcp)
fig9b<-mg_PlotMutiBoxplot(tcga.mcp[tcga.subtype$Samples,],
                          tcga.subtype$Cluster
                          , test_method = 'kruskal.test'
                          , add = 'boxplot'
                          , legend.pos = 'top'
                          , group_cols = ggsci::pal_lancet()(9)[c(1,4,3)]
                          , ylab = 'Enrichment Score',xangle = 30)
fig9b
write.table(tcga.mcp,'results/ana/09.subtype.immnu/tcga.mcp.txt',quote = F,sep='\t',row.names = T)

fig9ab<-mg_merge_plot(fig9a,fig9b,widths = c(1,2),nrow = 1,ncol = 2,labels = c('A','B'),align = 'hv')
fig9ab
#检查点
tcga.icg=immu_ICGs(tcga_tmp_log2_T)
fig9c=mg_PlotMutiBoxplot(tcga.icg[tcga.subtype$Samples,]
                         , group = tcga.subtype$Cluster
                         , legend.pos = 'right'
                         , add = 'boxplot'
                         , xangle = 30,
                         , ylab = 'log2(Gene Expression+1)'
                         , group_cols = ggsci::pal_lancet()(9)[c(1,4,3)]
                         , test_method = 'kruskal.test')
fig9c
write.table(tcga.icg,'results/ana/09.subtype.immnu/tcga.icg.txt',quote = F,row.names = T,sep='\t')
#TIDE
tcga_tide_dat <- t(scale(t(tcga_tmp_log2_T),scale = F))
write.table(tcga_tide_dat,
            file = 'results/ana/09.subtype.immnu/tcga_tide_dat.txt',
            quote = F, sep = '\t')

tcga_tide_dat<-read.csv('results/ana/09.subtype.immnu/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
tcga_tide_dat=tcga_tide_dat[rownames(tcga.subtype),]
tcga_tide_dat=cbind.data.frame(tcga_tide_dat,Cluster=tcga.subtype[rownames(tcga_tide_dat),]$Cluster)
head(tcga_tide_dat)
fig9d=list()

fig9d[[1]]=sig_boxplot(tcga_tide_dat[,c("Cluster","TIDE")],
                       leg = 'Cluster',ylab = 'TIDE',
                       palette =ggsci::pal_lancet()(9)[c(1,4,3)])
fig9d[[1]]

fig9d[[2]]=sig_boxplot(tcga_tide_dat[,c("Cluster","CD274")],
                       leg = 'Cluster',ylab = 'CD274',
                       palette =ggsci::pal_lancet()(9)[c(1,4,3)])
fig9d[[2]]
fig9d[[3]]=sig_boxplot(tcga_tide_dat[,c("Cluster","Dysfunction")],
                       leg = 'Cluster',ylab = 'Dysfunction',
                       palette =ggsci::pal_lancet()(9)[c(1,4,3)])
fig9d[[3]]

fig9d[[4]]=sig_boxplot(tcga_tide_dat[,c("Cluster","MDSC")],
                       leg = 'Cluster',ylab = 'MDSC',
                       palette =ggsci::pal_lancet()(9)[c(1,4,3)])
fig9d[[4]]
fig9d[[5]]<-plotMutiBar(table(tcga_tide_dat$Responder,tcga_tide_dat$Cluster),
                        ist = F,legTitle = 'Responder')
fig9d[[5]]
fig9d<-mg_merge_plot(fig9d,nrow = 1,ncol = 5,common.legend = T)
fig9<-mg_merge_plot(fig9ab,fig9c,fig9d,labels = c('','C','D'),nrow = 3,ncol = 1,heights = c(1.5,1,1.2))
ggsave('results/ana/09.subtype.immnu/Fig9.pdf',fig9,height = 12,width = 15)


##11、衰老相关风险模型的构建####
dir.create('results/ana/10.module')
#差异分析

gene_type=gene_type[gene_type$TYPE=='protein_coding',]
tcga.subtype$Cluster1=ifelse(tcga.subtype$Cluster=='clust1','clust1','clust23')
tcga.subtype$Cluster2=ifelse(tcga.subtype$Cluster=='clust2','clust2','clust13')
tcga.subtype$Cluster3=ifelse(tcga.subtype$Cluster=='clust3','clust3','clust12')

tcga.diff1<-mg_limma_DEG(exp = tcga_tmp_log2_T[as.character(unique(gene_type$SYMBOL)),tcga.subtype$Samples],
                         group = tcga.subtype$Cluster1,
                         ulab = 'clust1',dlab = 'clust23')
tcga.diff2<-mg_limma_DEG(exp = tcga_tmp_log2_T[as.character(unique(gene_type$SYMBOL)),tcga.subtype$Samples],
                         group = tcga.subtype$Cluster2,
                         ulab = 'clust2',dlab = 'clust13')
tcga.diff3<-mg_limma_DEG(exp = tcga_tmp_log2_T[as.character(unique(gene_type$SYMBOL)),tcga.subtype$Samples],
                         group = tcga.subtype$Cluster3,
                         ulab = 'clust3',dlab = 'clust12')
tcga.diff1<-tcga.diff1$DEG
tcga.diff2<-tcga.diff2$DEG
tcga.diff3<-tcga.diff3$DEG

fc_fit=1
p_fit=0.05
table(tcga.diff1$adj.P.Val<0.05)

tcga.diff1.fit<-tcga.diff1[which(abs(tcga.diff1$logFC)>fc_fit & tcga.diff1$P.Value<p_fit),]
table(tcga.diff1.fit$logFC>0)
# FALSE  TRUE 
#2   34
tcga.diff2.fit<-tcga.diff2[which(abs(tcga.diff2$logFC)>fc_fit & tcga.diff2$P.Value<p_fit),]
table(tcga.diff2.fit$logFC>0)
# FALSE  TRUE 
# 293   9 
tcga.diff3.fit<-tcga.diff3[which(abs(tcga.diff3$logFC)>fc_fit & tcga.diff3$P.Value<p_fit),]
table(tcga.diff3.fit$logFC>0)
# FALSE  TRUE 
# 26   218

write.table(tcga.diff1,'results/ana/10.module//tcga.diff1.txt',quote = F,sep='\t',row.names = T)
write.table(tcga.diff1,'results/files/tcga.diff1.txt',quote = F,sep='\t',row.names = T)

write.table(tcga.diff2,'results/ana/10.module/tcga.diff2.txt',quote = F,sep='\t',row.names = T)
write.table(tcga.diff2,'results/files/tcga.diff2.txt',quote = F,sep='\t',row.names = T)

write.table(tcga.diff3,'results/ana/10.module/tcga.diff3.txt',quote = F,sep='\t',row.names = T)
write.table(tcga.diff3,'results/files/tcga.diff3.txt',quote = F,sep='\t',row.names = T)

coxRun<-function(dat){
  library(survival)
  colnames(dat)=c('time','status','AS')  
  dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
  #print(nrow(dat))
  if(nrow(dat)<10){
    print(paste0('Sample Num is small:',nrow(dat)))
    return(c(NA,NA,NA,NA))
  }
  #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
  fmla <- as.formula("Surv(time, status) ~AS")
  if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
    cox <- survival::coxph(fmla, data = dat)
    re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
    return(re)
  }else{
    return(c(NA,NA,NA,NA))
  }
}
cox_batch<-function(dat,time,event){
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}

mg_risksocre.sig<-function(dat,os,os.time){
  #dat=tcga_cellage.exp_for
  #os=tcga_cli$OS
  #os.time=tcga_cli$OS.time
  tcga_dat_m<-cbind.data.frame(os=os,os.time=os.time,dat)
  tcga.cox <- cox_batch(t(tcga_dat_m[,3:ncol(tcga_dat_m)]),
                        time = tcga_dat_m$os.time,
                        event = tcga_dat_m$os)
  
  return(tcga.cox)
}
#lasso
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10,
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=sig.coef,lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
#多因素获取凤霞得分
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}

#单因素cox分析
diff.gene=unique(c(rownames(tcga.diff1.fit),rownames(tcga.diff2.fit),rownames(tcga.diff3.fit)))
length(diff.gene)#391
write.table(diff.gene,'results/files/all.diff.gene.txt',quote = F,row.names = F,sep='\t',col.names = F)

tcga_diff.exp=tcga_tmp_log2_T[diff.gene,tcga_subtype.cli$Samples]
rownames(tcga_diff.exp)=gsub('-','__',rownames(tcga_diff.exp))
tcga.cox <- mg_risksocre.sig(dat = t(tcga_diff.exp),
                             os = tcga_subtype.cli$OS,
                             os.time = tcga_subtype.cli$OS.time)
cox.p_fit=0.05
tcga.cox.fit<-tcga.cox[tcga.cox$p.value<cox.p_fit,]
dim(tcga.cox.fit)#69
tcga.cox$coef=log(tcga.cox$HR)
tcga.cox$Gene=rownames(tcga.cox)
tcga.cox$type=rep('None',nrow(tcga.cox))
tcga.cox$type[which(tcga.cox$p.value<cox.p_fit & tcga.cox$coef>0)]='Risk'
tcga.cox$type[which(tcga.cox$p.value<cox.p_fit & tcga.cox$coef<0)]='Protective'
table(tcga.cox$type)
fig10a <- ggplot(data = tcga.cox,
                aes(x = coef,
                    y = -log10(p.value)))+
  geom_point(alpha=0.4, size=3.5, aes(color=type))+
  scale_color_manual(values=c(mg_colors[2],'grey',mg_colors[1]),
                     limits = c("Protective",'None', "Risk"),name='State')+
  geom_hline(yintercept = -log10(cox.p_fit),lty=4,col="black",lwd=0.8)+
  ylab('-log10(pvalue)')+xlab('Cox coefficient')+
  theme_bw()+
  theme(axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
        legend.text=element_text(face="plain", family="Times", colour="black"),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black" ),#设置图例的总标题的字体属性
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill = NA, colour = NA)
  )
# +ggrepel::geom_text_repel(data=module.genes.cox_use[which(module.genes.cox_use$p.value<0.05),],aes(label=Gene))
fig10a

table(tcga.cox$type)
#  None Protective       Risk 
#322        6         63
dim(tcga.cox)
write.table(tcga.cox,'results/ana/10.module/tcga.cox.txt',quote = F,row.names = F,sep='\t')
write.table(tcga.cox,'results/files/tcga.cox.txt',quote = F,row.names = F,sep='\t')
options(ggrepel.max.overlaps=Inf)
tcga.lasso<-get_riskscore.lasso(dat = t(tcga_diff.exp[rownames(tcga.cox.fit),tcga_subtype.cli$Samples]),
                                os=tcga_subtype.cli$OS,
                                os.time = tcga_subtype.cli$OS.time,
                                labels=c('B','C'))

length(tcga.lasso$lasso.gene)#15
tcga.lasso$plot

ggsave('results/ana/10.module/Fig10.pdf',fig10a,height = 5,width = 5)
ggsave('results/ana/10.module/Fig10bc.pdf',tcga.lasso$plot,height = 5,width = 10)
tcga.lasso$lasso.gene
names(tcga.lasso$lasso.gene)

tcga.lasso$lambda.min
#  0.01727836
tcga.module.risk<-get_riskscore(dat=t(tcga_diff.exp[names(tcga.lasso$lasso.gene),]),
                                #dat=tcga_diff.exp[,rownames(tcga.cox.fit)],
                                os=tcga_subtype.cli$OS,
                                os.time=tcga_subtype.cli$OS.time,
                                step=F,
                                direction=c("both", "backward", "forward")[1])

length(tcga.module.risk$module.gene)#4
names(tcga.module.risk$module.gene)
#"MMP28" "CT83"  "ACP5"  "KRT6A"
tcga.module.risk$model
#0.094*MMP28+-0.346*CT83+0.483*ACP5+0.186*KRT6A

#柱状图
gene.coef=as.data.frame(tcga.module.risk$module.gene)
gene.coef=data.frame(Gene=rownames(gene.coef),Coef=gene.coef[,1])
gene.coef$Type=ifelse(gene.coef$Coef>0,'Risk','Protective')
gene.coef$Type=factor(gene.coef$Type,levels=c('Risk','Protective'))
library(dplyr)
fig10d=gene.coef %>% 
  ggplot(aes(reorder(Gene, Coef), Coef)) +
  geom_col(aes(fill = Type)) +
  coord_flip() +
  scale_fill_manual(values=ggsci::pal_lancet('lanonc',alpha =0.6)(9)[c(7,1)]) +
  coord_flip() +
  labs(x = "") +
  labs(y = "Cox coefficient") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1),legend.position="top")
fig10d
ggsave('results/ana/10.module/Fig10d.pdf',fig9d,height = 5,width = 5)


fig11a<-ggplotTimeROC(time = tcga.module.risk$result$time/365
                      ,status = tcga.module.risk$result$status
                      ,score = tcga.module.risk$result$riskscorez
                      ,mks = c(1,2,3,4,5))
fig11a
fig11b<-ggplotKMCox(dat = data.frame(tcga.module.risk$result$time/365,
                                     tcga.module.risk$result$status,
                                     #tcga.module.risk$result$Risk
                                     ifelse(tcga.module.risk$result$riskscorez>=0,'H','L')),
                    title = 'TCGA',add_text = '',
                    palette = ggsci::pal_aaas()(10)[c(2,1)],labs = c('High','Low'))
fig11b
tcga.module.risk$result$Risk=ifelse(tcga.module.risk$result$riskscorez>=0,'High','Low')
write.table(tcga.module.risk$result,'results/ana/10.module/tcga.module.risk.txt',quote = F,row.names = T,sep='\t')


#geo数据集验证
names(tcga.module.risk$module.gene)
gse58812.cli.module.risk<-get_riskscore(dat = t(GSE58812_exp[intersect(names(tcga.module.risk$module.gene),rownames(GSE58812_exp)),
                                                             GSE58812_cli$Samples]),
                                        os =GSE58812_cli$OS,
                                        os.time = GSE58812_cli$OS.time,
                                        step=F,direction=c("both", "backward", "forward")[1])

fig11c<-ggplotTimeROC(time = gse58812.cli.module.risk$result$time/365
                      ,status = gse58812.cli.module.risk$result$status
                      ,score = gse58812.cli.module.risk$result$riskscorez
                      ,mks = c(1,2,3,4,5))
fig11c
fig11d<-ggplotKMCox(dat = data.frame(gse58812.cli.module.risk$result$time/365,
                                     gse58812.cli.module.risk$result$status,
                                     ifelse(gse58812.cli.module.risk$result$riskscorez>=0,'H','L')),
                    #gse58812.cli.module.risk$result$Risk),
                    title = 'GSE58812',add_text = '',
                    palette = ggsci::pal_aaas()(10)[c(2,1)],labs = c('High','Low'))
fig11d

gse58812.cli.module.risk$result$Risk=ifelse(gse58812.cli.module.risk$result$riskscorez>=0,'High','Low')
write.table(gse58812.cli.module.risk$result,'results/ana/10.module/gse58812.cli.module.risk.txt',quote = F,sep='\t',row.names = T)
#cbio
cbio.cli.module.risk<-get_riskscore(dat = t(cbio_exp[intersect(names(tcga.module.risk$module.gene),rownames(cbio_exp)),
                                                     cbio_cli$Samples]),
                                        os =cbio_cli$OS,
                                        os.time = cbio_cli$OS.time,
                                        step=F,direction=c("both", "backward", "forward")[1])

fig11e<-ggplotTimeROC(time = cbio.cli.module.risk$result$time/365
                      ,status = cbio.cli.module.risk$result$status
                      ,score = cbio.cli.module.risk$result$riskscorez
                      ,mks = c(1,2,3,4,5))
fig11e
fig11f<-ggplotKMCox(dat = data.frame(cbio.cli.module.risk$result$time/365,
                                     cbio.cli.module.risk$result$status,
                                     ifelse(cbio.cli.module.risk$result$riskscorez>=0,'H','L')),
                    #cbio.cli.module.risk$result$Risk),
                    title = 'cbio',add_text = '',
                    palette = ggsci::pal_aaas()(10)[c(2,1)],labs = c('High','Low'))
fig11f

gse58812.cli.module.risk$result$Risk=ifelse(gse58812.cli.module.risk$result$riskscorez>=0,'High','Low')
write.table(gse58812.cli.module.risk$result,'results/ana/10.module/gse58812.cli.module.risk.txt',quote = F,sep='\t',row.names = T)

fig11<-mg_merge_plot(fig11a,fig11b,fig11c,fig11d,nrow = 2,ncol = 2,labels = c('A','B','C','D'))
fig11
ggsave('results/ana/10.module/Fig11.pdf',fig11,height = 9,width = 9)
#11、风险模型与临床表型的关系#######
dir.create('results/ana/11.module.cli')
head(tcga_subtype.cli)
tcga_subtype.cli.risk<-merge(tcga_subtype.cli,
                             tcga.module.risk$result[,c("Samples","riskscore","riskscorez","Risk")],
                             by='Samples')
head(tcga_subtype.cli.risk)
write.table(tcga_subtype.cli.risk,'results/ana/11.module.cli/tcga_subtype.cli.risk.txt',quote = F,row.names = F,sep='\t')
fig12a<-list()
fig12a[[1]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Cluster","riskscore")],
                         leg = 'Cluster',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[1]]
fig12a[[2]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("T.Stage","riskscore")],
                         leg = 'T.Stage',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[2]]
fig12a[[3]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("M.Stage","riskscore")],
                         leg = 'M.Stage',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[3]]
fig12a[[4]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("N.Stage","riskscore")],
                         leg = 'N.Stage',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[4]]
fig12a[[5]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Stage","riskscore")],
                         leg = 'Stage',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[5]]

fig12a[[6]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Age1","riskscore")],
                         leg = 'Age',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[6]]
fig12a[[7]]<-sig_boxplot(dat = tcga_subtype.cli.risk[,c("Event","riskscore")],
                         leg = 'Event',ylab = 'Senescore',
                         palette =ggsci::pal_aaas()(10))
fig12a[[7]]

fig12a1<-mg_merge_plot(fig12a,nrow = 2,ncol = 4)
fig12b<-list()
fig12b[[1]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$Cluster,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[1]]
fig12b[[2]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$T.Stage,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[2]]
fig12b[[3]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$M.Stage,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[3]]
fig12b[[4]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$N.Stage,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[4]]
fig12b[[5]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$Stage,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[5]]
fig12b[[6]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$Age1,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[6]]
fig12b[[7]]<-plotMutiBar(dat = table(tcga_subtype.cli.risk$Event,tcga_subtype.cli.risk$Risk),
                         ist = T,legTitle = 'Risk')
fig12b[[7]]
fig12b1<-mg_merge_plot(fig12b,nrow = 2,ncol = 4)

fig12<-mg_merge_plot(fig12a1,fig12b1,nrow = 2,ncol = 1,labels = c('A','B'),heights = c(1,1))
ggsave(plot = fig12,filename = 'results/ana/11.module.cli/Fig12.pdf',height = 15,width = 15)



#12单因素和多因素+列线图#####
dir.create('results/ana/12.module.sigmut')
library(forestplot)
library(survcomp)
tcga_cox_datas=tcga_subtype.cli.risk
table(tcga_cox_datas$Cluster)
table(tcga_cox_datas$Gender)
table(tcga_cox_datas$T.Stage)
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage=='T1'|tcga_cox_datas$T.Stage=='T2']<-'T1+T2'
tcga_cox_datas$T.Stage[tcga_cox_datas$T.Stage=='T3'|tcga_cox_datas$T.Stage=='T4']<-'T3+T4'

table(tcga_cox_datas$M.Stage)

table(tcga_cox_datas$N.Stage)
tcga_cox_datas$N.Stage[tcga_cox_datas$N.Stage=='N1'|tcga_cox_datas$N.Stage=='N2'|tcga_cox_datas$N.Stage=='N3']<-'N1+N2+N3'

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'


tcga_cox_datas$RiskType=tcga_cox_datas$Risk
table(tcga_cox_datas$RiskType)
tcga_cox_datas$RiskType <- ifelse(tcga_cox_datas$RiskType == 'High', '1High', '0Low')

write.table(tcga_cox_datas,'results/ana/12.module.sigmut/tcga_cox_datas.txt',quote = F,sep = '\t',row.names = F)
# 单因素分析结果 
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat


A3_T_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.Stage,
                              data=tcga_cox_datas))
A3_T_sig_cox_dat <- data.frame(Names=rownames(A3_T_sig_cox[[8]]),
                               HR = round(A3_T_sig_cox[[7]][,2],3),
                               lower.95 = round(A3_T_sig_cox[[8]][,3],3),
                               upper.95 = round(A3_T_sig_cox[[8]][,4],3),
                               p.value=round(A3_T_sig_cox[[7]][,5],3))
A3_T_sig_cox_dat

A4_N_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.Stage,
                              data=tcga_cox_datas))
A4_N_sig_cox_dat <- data.frame(Names=rownames(A4_N_sig_cox[[8]]),
                               HR = round(A4_N_sig_cox[[7]][,2],3),
                               lower.95 = round(A4_N_sig_cox[[8]][,3],3),
                               upper.95 = round(A4_N_sig_cox[[8]][,4],3),
                               p.value=round(A4_N_sig_cox[[7]][,5],3))
A4_N_sig_cox_dat


Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


RiskType_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~riskscorez,
                                  data=tcga_cox_datas))
RiskType_sig_cox_dat <- data.frame(Names=rownames(RiskType_sig_cox[[8]]),
                                   HR = round(RiskType_sig_cox[[7]][,2],3),
                                   lower.95 = round(RiskType_sig_cox[[8]][,3],3),
                                   upper.95 = round(RiskType_sig_cox[[8]][,4],3),
                                   p.value=round(RiskType_sig_cox[[7]][,5],3))
RiskType_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     A3_T_sig_cox_dat,
                     A4_N_sig_cox_dat,
                     Stage_sig_cox_dat,
                     RiskType_sig_cox_dat)

data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "T.Stage",
                        "N.Stage",
                        "Stage",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)

pdf('results/ana/12.module.sigmut/Fig13A.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 8,lineheight = 10)
dev.off()
#多因素分析结果
#muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender+T.Stage+N.Stage+Grade+Stage+RiskType, data=tcga_cox_datas))
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+T.Stage+N.Stage+Stage+riskscorez, data=tcga_cox_datas))

muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
rownames(data.muti) <- c("Age",
                         "T.Stage",
                         "N.Stage",
                         "Stage",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)

pdf('results/ana/12.module.sigmut/Fig13B.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10)
dev.off()
#列线图和DCA
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

pdf('results/ana/12.module.sigmut/fig13C.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$riskscorez,
                                Stage=tcga_cox_datas$Stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))

#AUC
library("timeROC")
tcga_cox_auc=tcga_cox_datas
tcga_cox_auc$T.Stage=as.numeric(as.factor(tcga_cox_auc$T.Stage))
tcga_cox_auc$N.Stage=as.numeric(as.factor(tcga_cox_auc$N.Stage))
tcga_cox_auc$M.Stage=as.numeric(as.factor(tcga_cox_auc$M.Stage))

tcga_cox_auc$Stage=as.numeric(as.factor(tcga_cox_auc$Stage))
tcga_cox_auc$riskscore=as.numeric(tcga_cox_auc$riskscore)
tcga_cox_auc$Age=as.numeric(tcga_cox_auc$Age)


ROC.DSST.Age=timeROC(T=tcga_cox_auc$OS.time/365,
                     delta=tcga_cox_auc$OS,
                     marker=tcga_cox_auc$Age,
                     cause=1,weighting="marginal",
                     times=c(1,2,3,4,5),
                     iid=TRUE)
ROC.DSST.T.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$T.Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.M.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$M.Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)
ROC.DSST.N.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                         delta=tcga_cox_auc$OS,
                         marker=tcga_cox_auc$N.Stage,
                         cause=1,weighting="marginal",
                         times=c(1,2,3,4,5),
                         iid=TRUE)

ROC.DSST.Stage=timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$Stage,
                       cause=1,weighting="marginal",
                       times=c(1,2,3,4,5),
                       iid=TRUE)

ROC.DSST.Risk=timeROC(T=tcga_cox_auc$OS.time/365,
                      delta=tcga_cox_auc$OS,
                      marker=tcga_cox_auc$riskscorez,
                      cause=1,weighting="marginal",
                      times=c(1,2,3,4,5),
                      iid=TRUE)

ROC.DSST.Nomo<-timeROC(T=tcga_cox_auc$OS.time/365,
                       delta=tcga_cox_auc$OS,
                       marker=tcga_cox_auc$riskscore,
                       other_markers=as.matrix(tcga_cox_auc[,c("Stage")]),
                       cause=1,
                       weighting="cox",
                       times=c(1,2,3,4,5),
                       iid=F)
ROC.DSST.Age$AUC
ROC.DSST.T.Stage$AUC
ROC.DSST.M.Stage$AUC
ROC.DSST.N.Stage$AUC
ROC.DSST.Stage$AUC
ROC.DSST.Risk$AUC
ROC.DSST.Nomo$AUC

scales::show_col(mg_colors[c(1,10:12,4,5,7:9)])

pdf('results/ana/12.module.sigmut/Fig13J.pdf',height = 5,width = 6)
plotAUCcurve(ROC.DSST.Nomo,conf.int=F,col=mg_colors[1])
plotAUCcurve(ROC.DSST.Age,conf.int=F,col=mg_colors[2],add=TRUE)
plotAUCcurve(ROC.DSST.T.Stage,conf.int=F,col=mg_colors[3],add=TRUE)
plotAUCcurve(ROC.DSST.M.Stage,conf.int=F,col=mg_colors[4],add=TRUE)
plotAUCcurve(ROC.DSST.N.Stage,conf.int=F,col=mg_colors[5],add=TRUE)
plotAUCcurve(ROC.DSST.Stage,conf.int=F,col=mg_colors[6],add=TRUE)
plotAUCcurve(ROC.DSST.Risk,conf.int=F,col=mg_colors[7],add=TRUE)
legend("topright",c("Nomogram","Age","T.Stage",'M.Stage','N.Stage','Stage','RiskScore')
       ,col=mg_colors[c(1:7)],lty=1,lwd=2)

dev.off()

#13、风险模型的免疫特征####
dir.create('results/ana/13.module.immnu')
head(tcga.mcp)
ggsci::pal_aaas()(10)[c(2,1)]
fig14a<-mg_PlotMutiBoxplot(data = tcga.mcp[tcga_subtype.cli.risk$Samples,],
                           group = tcga_subtype.cli.risk$Risk,
                           group_cols = ggsci::pal_aaas()(10)[c(2,1)],
                           test_method = 'wilcox.test',
                           xangle = 30,
                           add = 'boxplot', 
                           legend.pos = 'top',ylab = 'Immnune score')
fig14a
write.table(cbind.data.frame(tcga.mcp[tcga_subtype.cli.risk$Samples,],Risk=tcga_subtype.cli.risk$Risk),'results/ana/13.module.immnu/tcga.mcp.risk.txt',quote = F,row.names = T,sep='\t')
#Estimate
fig14b<-mg_PlotMutiBoxplot(data = tcga_esti[tcga_subtype.cli.risk$Samples,],
                           group = tcga_subtype.cli.risk$Risk,
                           group_cols = ggsci::pal_aaas()(10)[2:1],
                           test_method = 'wilcox.test',
                           xangle = 30,
                           add = 'boxplot', 
                           legend.pos = 'top',ylab = 'Immnune score')
fig14b
write.table(cbind.data.frame( tcga_esti[tcga_subtype.cli.risk$Samples,],Risk=tcga_subtype.cli.risk$Risk),'results/ana/13.module.immnu/tcga.esti.risk.txt',quote = F,row.names = T,sep='\t')
fig14ab<-mg_merge_plot(fig14a,fig14b,labels = c('A','B'),widths = c(3,1),nrow = 1,ncol = 2)

#高低风险组与免疫检查点
fig14c<-mg_PlotMutiBoxplot(data = tcga.icg[tcga_subtype.cli.risk$Samples,],
                           group = tcga_subtype.cli.risk$Risk,
                           group_cols = ggsci::pal_aaas()(10)[2:1],
                           test_method = 'wilcox.test',
                           xangle = 30,
                           add = 'boxplot', 
                           legend.pos = 'right',ylab = 'log2(Gene Expression +1)')
fig14c
fig14abc<-mg_merge_plot(fig14ab,fig14c,nrow = 2,ncol = 1,labels = c('','C'),heights = c(1,0.8))
ggsave('results/ana/13.module.immnu/fig14abc.pdf',fig14abc,height = 9,width = 15)
#相关性分析
cib.risk.cor<-data.frame(tcga.mcp[tcga_subtype.cli.risk$Samples,],RiskScore=tcga_subtype.cli.risk$riskscore)
cib.risk.cor.res <- Hmisc::rcorr(as.matrix(cib.risk.cor), 
                                 type = 'pearson')
cib.risk.cor.res$r
cib.risk.cor.res$P[is.na(cib.risk.cor.res$P)] <- 0
library(corrplot)
col1 = colorRampPalette(c('blue', 'white', 'red'))
pdf('results/ana/13.module.immnu/Fig14D.pdf',height = 7,width = 7)
corrplot(as.matrix(cib.risk.cor.res$r), 
         p.mat = as.matrix(cib.risk.cor.res$P),
         mar = c(0,0,1,1),
         col=col1(100),
         tl.srt = 90,
         tl.cex = 1,
         tl.col = 'black',
         tl.offset = 0.5,
         cl.pos = c("b","r","n")[1], 
         #cl.lim=c(-0.3,0.3),
         cl.align.text = 'l',
         cl.length = 5,
         cl.ratio = 0.1,
         cl.cex = 0.8,
         addgrid.col = 'white',
         method = 'color',
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,
         is.corr=T,
         xpd=T
)
dev.off()
tcga.risk.immnu<-data.frame(tcga_esti[tcga_subtype.cli.risk$Samples,],RiskScore=tcga_subtype.cli.risk$riskscore)
head(tcga.risk.immnu)


fig14E1<-ggplot(tcga.risk.immnu, aes(x = RiskScore, y = ImmuneScore))+
  geom_point()+
  geom_smooth(method = "lm", color = "blue", fill = "lightgray")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+xlab('RiskScore')+ylab('ImmuneScore')
fig14E1
fig14E2<-ggplot(tcga.risk.immnu, aes(x = RiskScore, y = StromalScore))+
  geom_point()+
  geom_smooth(method = "lm", color = "blue", fill = "lightgray")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+xlab('RiskScore')+ylab('StromalScore')
fig14E2
fig14E3<-ggplot(tcga.risk.immnu, aes(x = RiskScore, y = ESTIMATEScore))+
  geom_point()+
  geom_smooth(method = "lm", color = "blue", fill = "lightgray")+
  stat_cor(method = "pearson")+
  theme_bw()+
  theme(panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+xlab('RiskScore')+ylab('ESTIMATEScore')
fig14E3
fig14E<-mg_merge_plot(fig14E1,fig14E2,fig14E3,nrow =1,ncol = 3)
ggsave('results/ana/13.module.immnu/fig14E.pdf',fig14E,height = 4,width = 12)
#IC50药物
#IC50药物
# library(pRRophetic)
# library(ggplot2)
# ## Cisplatin,顺铂
# set.seed(12345)
# predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_tmp_log2_T[,tcga.subtype$Samples])
#                                               , "Cisplatin"
#                                               , "urogenital_system"
#                                               , selection=1
#                                               ,dataset = "cgp2016")
# predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)
# 
# tcga_durg_ic50_res <- predictedPtype_Cisplatin
# 
drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
# dim(tcga_durg_ic50_res)
# for (drug in drugs) {
#   print(drug)
#   set.seed(12345)
#   tmpic50 <- pRRopheticPredict(as.matrix(tcga_tmp_log2_T[,tcga.subtype$Samples])
#                                , drug
#                                , "urogenital_system"
#                                , selection=1
#                                , dataset = "cgp2016")
#   tmpic50 <- data.frame(tmpic50)
#   colnames(tmpic50) <- drug
#   tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
# }
# 
# dim(tcga_durg_ic50_res)
# colnames(tcga_durg_ic50_res)
# colnames(tcga_durg_ic50_res) <- gsub('predictedPtype_', '', colnames(tcga_durg_ic50_res))
# save(tcga_durg_ic50_res,file = 'tcga_durg_ic50_res.RData')
load('tcga_durg_ic50_res.RData')
dim(tcga_durg_ic50_res)
tcga_durg_ic50_res2=cbind.data.frame(tcga_durg_ic50_res[tcga_subtype.cli.risk$Samples,],
                                     Risk=tcga_subtype.cli.risk$Risk,
                                     RiskScore=tcga_subtype.cli.risk$riskscore)
head(tcga_durg_ic50_res2)
write.table(tcga_durg_ic50_res2,'results/ana/13.module.immnu/tcga_durg_ic50_res2.txt',quote = F,row.names = T,sep='\t')

fig14f<-list()

fig14f[[1]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[6])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[6],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[1]]

fig14f[[2]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[10])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[10],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[2]]
#9，18，19，21，28，29,33,35,48

fig14f[[3]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[22])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[22],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[3]]

fig14f[[4]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[36])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[36],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[4]]

fig14f[[5]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[46])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[46],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[5]]
fig14f[[6]]<-sig_boxplot(dat = tcga_durg_ic50_res2[,c('Risk',drugs[9])],
                         leg = 'Risk',ylab = paste0('Estimated IC50 (',drugs[9],')'),
                         palette = ggsci::pal_aaas()(10)[2:1])
fig14f[[6]]

fig14f<-mg_merge_plot(fig14f,nrow = 1,ncol = 6,common.legend = T)
ggsave('results/ana/13.module.immnu/fig14f.pdf',fig14f,height = 3,width = 15)
#TIDE
# #Tide的差异分析
# head(tcga_tide_dat)
# tcga_tide_dat.risk=cbind.data.frame(tcga_tide_dat[tcga_subtype.cli.risk$Samples,],
#                                     Risk=tcga_subtype.cli.risk$Risk,
#                                     RiskScore=tcga_subtype.cli.risk$riskscore)
# head(tcga_tide_dat.risk)
# write.table(tcga_tide_dat.risk,'results/ana/13.module.immnu/tcga_tide_dat.risk.txt',quote = F,row.names = T,sep='\t')
# 
# #tide.selected=c('MDSC','IFNG','TAM.M2','Exclusion','Dysfunction','TIDE')
# fig14g=list()
# fig14g[[1]]=sig_boxplot(tcga_tide_dat.risk[,c("Risk","TIDE")],
#                     leg = 'Risk',ylab = 'TIDE',
#                     palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[1]]
# fig14g[[2]]=sig_boxplot(tcga_tide_dat.risk[,c("Risk","IFNG")],
#                     leg = 'Risk',ylab = 'IFNG',
#                     palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[2]]
# fig14g[[3]]=sig_boxplot(tcga_tide_dat.risk[,c("Risk","MDSC")],
#                     leg = 'Risk',ylab = 'MDSC',
#                     palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[3]]
# 
# fig14g[[4]]=sig_boxplot(tcga_tide_dat.risk[,c("Risk","Exclusion")],
#                     leg = 'Risk',ylab = 'Exclusion',
#                     palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[4]]
# fig14g[[5]]<-sig_boxplot(tcga_tide_dat.risk[,c("Risk","Dysfunction")],
#                      leg = 'Risk',ylab = 'Dysfunction',
#                      palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[5]]
# fig14g[[6]]<-sig_boxplot(tcga_tide_dat.risk[,c("Risk","TAM.M2")],
#                      leg = 'Risk',ylab = 'TAM.M2',
#                      palette =ggsci::pal_aaas()(10)[2:1] )
# fig14g[[6]]
# fig14g<-mg_merge_plot(fig14g,nrow = 1,ncol = 6,common.legend = T)
# fig14g



save.image('Project_002.RData')

