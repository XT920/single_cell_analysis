# single_cell_analysis
10X Genomics single-cell data processing and analysis flow based on Seurat

# 单细胞RNA测序标准分析流程

基于Seurat的10X Genomics单细胞数据处理与分析流程

## 项目概述

本项目提供了一个完整的单细胞RNA测序数据分析流程，专门针对10X Genomics平台生成的数据。使用R语言和Seurat包实现从原始数据到细胞注释的完整分析，包括质量控制、标准化、聚类、降维和细胞类型识别等关键步骤。

## 功能特性

-  完整的质量控制与数据过滤
-  高变基因识别与数据标准化
- PCA线性降维与可视化
- UMAP/t-SNE非线性降维
- 细胞聚类与差异表达分析
- 基于SingleR的自动细胞类型注释
- 结果保存与可视化输出

## 安装要求

### 系统要求

- R (版本 ≥ 4.0.0)
- 至少8GB内存（推荐16GB以上处理大型数据集）

### R包依赖

```R
install.packages(c("Seurat", "patchwork", "dplyr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleR", "celldex"))
```

## 项目结构

```bash
single_cell_pipeline/
├── data/
│   └── filtered_feature_bc_matrix/  # 10X Genomics数据文件
│       ├── barcodes.tsv.gz
│       ├── features.tsv.gz
│       └── matrix.mtx.gz
├── scripts/
│   └── single_cell_analysis.R      # 主分析脚本
├── results/                        # 分析结果输出
├── figs/                           # 生成的可视化图表
└── README.md
```

## 使用方法

### 1. 准备数据

将10X Genomics输出文件放置在`data/filtered_feature_bc_matrix/`目录下，确保包含以下三个文件：

- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

### 2. 运行分析

在RStudio或R命令行中执行主脚本：

```R
source("scripts/single_cell_analysis.R")
```

### 3. 调整参数

根据您的数据特性，可能需要调整以下参数：

- 线粒体基因模式（小鼠数据使用`"^mt-"`）
- QC阈值（`nFeature_RNA`和`percent.mt`）
- PCA分析维度数
- 聚类分辨率参数

### 4. 查看结果

分析完成后，结果将保存在：

- `results/`：包含差异表达分析结果等数据
- `figs/`：包含所有生成的可视化图表
- `pbmc10k_final.rds`：完整的Seurat对象，可用于后续分析

## 分析流程步骤

1. **数据加载**：读取10X Genomics格式数据
2. **质量控制**：计算质控指标并过滤低质量细胞
3. **数据标准化**：使用LogNormalize方法标准化表达数据
4. **特征选择**：识别高变基因用于下游分析
5. **数据缩放**：缩放数据使基因表达具有可比性
6. **线性降维**：执行PCA分析降维
7. **细胞聚类**：基于图论的细胞聚类方法
8. **非线性降维**：UMAP/t-SNE可视化
9. **差异表达分析**：识别簇特异性标记基因
10. **细胞类型注释**：使用SingleR进行自动细胞类型鉴定
11. **结果保存**：保存分析结果和可视化图表

## 自定义分析

用户可以根据需要修改脚本中的以下部分：

### 调整QC阈值

```R
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                        nFeature_RNA < 2500 & 
                        percent.mt < 5)  # 修改这些阈值
```

### 更改聚类分辨率

```R
pbmc <- FindClusters(pbmc, resolution = 0.5)  # 增加分辨率可获得更多簇
```

### 使用不同的参考数据集进行细胞注释

```R
# 可使用不同的参考数据集
# hpca.se <- HumanPrimaryCellAtlasData()
# blueprint <- BlueprintEncodeData()
# dice <- DatabaseImmuneCellExpressionData()
```

## 常见问题解答

### Q: 处理小鼠数据需要注意什么？

A: 需要将线粒体基因模式从`"^MT-"`改为`"^mt-"`：

```R
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
```

### Q: 如何增加聚类数量？

A: 提高FindClusters函数中的resolution参数：

```R
pbmc <- FindClusters(pbmc, resolution = 0.8)  # 更高的分辨率产生更多簇
```

### Q: 分析需要多长时间？

A: 处理时间取决于数据规模。1万个细胞通常需要15-30分钟，10万个细胞可能需要几小时。

## 结果解释

分析完成后，您可以查看以下关键结果：

1. **质量控制图**：评估细胞质量和过滤效果
2. **PCA和UMAP图**：可视化细胞在降维空间中的分布
3. **差异表达表**：识别每个细胞簇的特异性标记基因
4. **细胞类型注释**：基于参考数据集的自动细胞类型识别
5. **特征表达图**：可视化特定基因在不同细胞簇中的表达

## 参考文献

- Butler, A., Hoffman, P., Smibert, P. et al. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nature Biotechnology.
- Stuart, T., Butler, A., Hoffman, P. et al. (2019). Comprehensive Integration of Single-Cell Data. Cell.
- Aran, D., Looney, A.P., Liu, L. et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology.
