# 设置项目目录
project_path <- "D:/Dev/projects/tgga"
if(!dir.exists(project_path)){
  stop("工作目录不存在")
}
setwd(project_path)

# 设置项目
item <- "chol"
######表达矩阵批量读入与合并######
# 先读入两个tsv检查gene_id和gene_name在不同样本中是否仍是一致的:
tsv1 <- "/0c9aba87-406f-4789-9275-e1d25bb3aea7/72fffe3e-d4fb-4862-ba01-e91289ed94ef.rna_seq.augmented_star_gene_counts.tsv"
tsv2 <- "/0e0dcae8-4cc2-4ba5-bc0d-07c9dbd0e5a2/d94cd170-9dff-4b44-a7a8-6d1c277b0f8b.rna_seq.augmented_star_gene_counts.tsv"
file1 <- paste(project_path,"/datasets/",item,tsv1,sep="")
file2 <- paste(project_path,"/datasets/",item,tsv2,sep="")
if(!file.exists(file1)){
  stop("file1不存在")
}
if(!file.exists(file2)){
  stop("file2不存在")
}
table1 <- data.table::fread(file1)
table2 <- data.table::fread(file2)
# 返回逻辑值TRUE，确认一致
gene_id_res <- identical(table1$gene_id,table2$gene_id)
if(!gene_id_res){
  stop("file1和file2的gene_id不一致")
}
# 返回逻辑值TRUE，确认一致
gene_name_res <- identical(table1$gene_name,table2$gene_name)
if(!gene_name_res){
  stop("file1和file2的gene_name不一致")
}

## 获取gene_id列
## 注意，data.table::fread方式读入文件并不能指定行名，所以才使用这种方法
x1 <- data.table::fread(file1,select = c("gene_id"))
head(x1)