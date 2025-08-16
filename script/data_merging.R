# 设置项目目录
project_path <- "D:/Dev/projects/tgga"
if(!dir.exists(project_path)){
  stop("工作目录不存在")
}
setwd(project_path)

# 设置项目
item <- "chol"
item_path <- paste(project_path,"/datasets/",item,sep="")

json_path <- paste(project_path,'/json/',item,sep="")
######表达矩阵批量读入与合并######
# 先读入两个tsv检查gene_id和gene_name在不同样本中是否仍是一致的:
tsv1 <- "/0c9aba87-406f-4789-9275-e1d25bb3aea7/72fffe3e-d4fb-4862-ba01-e91289ed94ef.rna_seq.augmented_star_gene_counts.tsv"
tsv2 <- "/0e0dcae8-4cc2-4ba5-bc0d-07c9dbd0e5a2/d94cd170-9dff-4b44-a7a8-6d1c277b0f8b.rna_seq.augmented_star_gene_counts.tsv"
file1 <- paste(item_path,tsv1,sep="")
file2 <- paste(item_path,tsv2,sep="")
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

# 批量读取所有tsv格式后缀文件的文件名：
# *表示任意前缀，$表示固定后缀
file_list <- dir(item_path,
               pattern = "*.rna_seq.augmented_star_gene_counts.tsv$",
               recursive = T)
head(file_list)

# 创建自定义函数,用于批量读入所有tsv文件的unstranded列（counts）：
exp_dt <- function(x){
  result <- data.table::fread(file.path(item_path,x),
                              select = c("unstranded"))
  return(result)
}
# 读取所有.tsv的unstranded列（将file_list转换为list，并对list中的每一个元素都应用函数exp_dt）
exp <- lapply(file_list,exp_dt)
# 将exp按列合并，并将list转化为data.table
exp <- do.call(cbind,exp)
exp[1:6,1:6]

# 添加行名(gene_id),所有表达矩阵的gene_id是相同的：
# data.table格式不能自定行名，因此我们先转换为数据框
exp <- as.data.frame(exp)
rownames(exp) <- x1$gene_id
exp[1:8,1:4]

# 列名替换为样本名
# 读取json文件（读取函数fromJSON来自R自带R包jsonlite）：
json_file <- paste(json_path,"/metadata.cart.2025-08-13.json",sep="")
if(!file.exists(json_file)){
  stop("JSON文件不存在")
}
meta_dt <- jsonlite::fromJSON(json_file)
# 查看这个list中第三列中的第一个list(包含数据框):
meta_dt[[3]][[1]]

id <- meta_dt$associated_entities
id[[1]]
# 这是我们需要的样本
id[[1]][,1]
# 这是我们需要用来匹配的文件名
head(meta_dt[[4]])
