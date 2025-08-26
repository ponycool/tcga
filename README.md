# TCGA数据处理工具

## 项目简介
本项目旨在简化TCGA(The Cancer Genome Atlas)数据的下载与基因表达矩阵的合成流程，为生物信息学研究人员提供高效、可重复的数据预处理解决方案。通过自动化脚本实现从数据下载到矩阵合成的完整工作流，减少手动操作错误，提高研究效率。

## 目录结构
```
TCGA/
├── GDC Data Transfer Tool/  # 存放GDC官方数据传输工具
├── data/                    # 保存处理后的基因表达矩阵等结果数据
├── datasets/                # 存储从TCGA下载的原始样本数据
├── json/                    # 保存数据下载相关的JSON元数据文件
├── script/                  # 包含数据处理的R语言脚本
└── README.md                # 项目说明文档
```

## 安装说明

### 前提条件
- R (>= 4.0)
- GDC Data Transfer Tool

### 项目获取
```bash
git clone https://github.com/ponycool/tcga.git
cd tcga
```

### 依赖安装
在R环境中安装必要的依赖包：
```r
install.packages(c("stringr", "jsonlite", "dplyr", "SummarizedExperiment"))
```

## 使用指南

### 1. 数据下载
#### 准备manifest文件
1. 登录[GDC Data Portal](https://portal.gdc.cancer.gov/)
2. 选择感兴趣的数据集，生成manifest文件
3. 将manifest文件保存到`json/`目录

#### 执行数据下载
```bash
# 使用GDC Data Transfer Tool下载数据
./GDC\ Data\ Transfer\ Tool/gdc-client download -m json/manifest.txt -d datasets/
```

### 2. 基因表达矩阵合成
运行R脚本处理下载的数据并生成基因表达矩阵：
```bash
Rscript script/data_mergin.R
```

## 脚本说明
- `data_mergin.R`: 实现数据解析、标准化和矩阵合成功能

## 贡献指南
1. Fork本仓库
2. 创建特性分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add some amazing feature'`)
4. 推送到分支 (`git push origin feature/amazing-feature`)
5. 打开Pull Request

## 许可证
本项目采用MIT许可证 - 详见[LICENSE](LICENSE)文件

## 致谢
- [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)项目提供数据支持