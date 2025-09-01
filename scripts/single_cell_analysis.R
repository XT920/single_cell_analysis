# 单细胞RNA-seq标准分析流程 (10X Genomics数据)
# 使用Seurat包进行分析
# 作者：许童


# 1. 安装和加载所需包 ----------------------------------------------------
install.packages('Seurat')
install.packages('patchwork')
install.packages('dplyr')
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleR")
BiocManager::install("celldex")

library(Seurat)
library(patchwork)
library(dplyr)
library(SingleR)
library(celldex)

# 2. 加载10X Genomics数据 -------------------------------------------------
pbmc.data <- Read10X(data.dir = "filtered_feature_bc_matrix/")

# 创建Seurat对象
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc10k", 
                           min.cells = 3, 
                           min.features = 200)

# 3. 质量控制 -------------------------------------------------------------
# 计算线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 可视化QC指标
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 过滤细胞
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)

# 4. 数据标准化 -----------------------------------------------------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 5. 特征选择 ------------------------------------------------------------
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 识别前10个高变基因
top10 <- head(VariableFeatures(pbmc), 10)

# 可视化高变基因
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 6. 数据缩放 ------------------------------------------------------------
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

# 7. 线性降维(PCA) -------------------------------------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# 可视化PCA结果
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")

# 8. 确定数据集维度 ------------------------------------------------------
# 通过JackStraw和ElbowPlot确定显著的主成分
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

# 9. 细胞聚类 ------------------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 10. 非线性降维(UMAP/t-SNE) ---------------------------------------------
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# 11. 寻找差异表达特征 ---------------------------------------------------
# 识别cluster 0的标记基因
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)

# 识别所有cluster的标记基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# 12. 细胞类型注释 -------------------------------------------------------
# 使用SingleR进行自动注释
hpca.se <- HumanPrimaryCellAtlasData()
pbmc.sce <- as.SingleCellExperiment(pbmc)
annotations <- SingleR(test = pbmc.sce, ref = hpca.se, labels = hpca.se$label.main)

# 将注释结果添加到Seurat对象
pbmc[["celltype"]] <- annotations$labels[match(rownames(pbmc@meta.data), rownames(annotations))]

# 可视化注释结果
DimPlot(pbmc, group.by = "celltype", label = TRUE)

# 13. 保存结果 -----------------------------------------------------------
saveRDS(pbmc, file = "pbmc10k_final.rds")

# 14. 高级可视化 ---------------------------------------------------------
# 特征表达可视化
FeaturePlot(pbmc, 
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A"), 
            reduction = "umap")

# 热图展示顶级标记基因
top10_markers <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10_markers$gene)
