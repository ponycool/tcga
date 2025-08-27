# 设置项目目录
project_path <- "D:/Dev/projects/tcga"
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

# 读取所有样本id，返回一个包含所有id名的向量（和前面lapply的区别，lapply返回的是list）
ids <- sapply(id,function(x){x[,1]}) 
# 生成一个包含有样本id和文件名这两列对应关系的数据框，文件名即为表达矩阵的文件名
ids2 <- data.frame(file_name = meta_dt$file_name,
                   id = ids) 
head(ids2)

# 观察文件名结构：
ids2$file_name[1:2]
head(file_list)[1:2]

# 用"/"做一个字符串拆分：把表达矩阵的文件名(file_list)的后半段和file_name做匹配：
library(stringr)
# 以/做拆分，simplify = T把结果返回成矩阵
file_list_dt <- str_split(file_list,"/",simplify = T)
head(file_list_dt)

# 取矩阵第二列：
file_list_dt <- file_list_dt[,2]
# 用百分百in检测是否拆分后的文件名都完全包含在file_name中(此时顺序不同，因此不能用==)
# 全部返回TRUE，则拆分无误
res <- file_list_dt %in% ids2$file_name
if(all(res)){
  message("SUCCESS 文件名全部包含在ids的filename列中")
}else{
  stop("FAILED 文件名未全部包含在ids的filename列中")
}

# 用match函数调整file_name的顺序和file_list_dt一致（因为我们合并的表达矩阵是按照file_list_dt的顺序进行合并的）
ids3 <- ids2[match(file_list_dt,ids2$file_name),]
# 查看并检查顺序：
head(ids3$file_name)
head(file_list_dt)
res <- identical(ids3$file_name,file_list_dt)
# 返回TRUE，确认无误
if(isTRUE(res)){
  message('SUCCESS 检查通过')
}else{
  stop("FAILED 检查未通过")
}

# 修改列名为样本名：
colnames(exp) <- ids3$id
exp[1:8,1:2]
message("表达矩阵合并完成！")

## 随便选一个tsv取gene_name这一列（这里记得转换为数据框，不然后面合并会丢失exp的行名）：
list <- as.data.frame(data.table::fread(file1,
                                        select = c("gene_name")))
exp2 <- cbind(list,exp)
#包含gene_name列的矩阵完成
exp2[1:10,1:3]


################# 方法一：
# 不需要区分九种不同的lncRNA了
lncRNA <- c("lncRNA")
mRNA <- c("protein_coding")

list2 <- as.data.frame(data.table::fread(file1,
                                         select = c("gene_id","gene_name","gene_type")))
lncRNA_list <- list2[list2$gene_type %in% lncRNA,]
mRNA_list <- list2[list2$gene_type %in% mRNA,]

# 现16901个，更新前14826个
length(lncRNA_list$gene_id)
# 现19962个，更新前19814个
length(mRNA_list$gene_id)

# v33版的GENECODE比v22版鉴定到的mRNA和lncRNA数量都更多了。

# 分别取表达矩阵(exp2)和mRNA/lncRNA_list的gene_id交集：
# mRNA交集：
mRNA_int <- intersect(
  rownames(exp2),
  mRNA_list$gene_id
)
# lncRNA交集：
lncRNA_int <- intersect(
  rownames(exp2),
  lncRNA_list$gene_id
)

mRNA_exp <- exp2[mRNA_int,]
lncRNA_exp <- exp2[lncRNA_int,]
###### 方法一结束

############# 方法二(直接用%in%判断并筛选即可)：
# 这个方法就是我们先把表达矩阵和gene-name、gene_type三者合并，通过gene_type筛选，最后把不需要的列去掉。
list22 <- as.data.frame(data.table::fread(file1,
                                          select = c("gene_name","gene_type")))
# 把gene_type列合并进来
exp3 <- cbind(list22,exp)
exp3[1:8,1:3]

# 下面直接筛选即可：
mRNA_exp2 <- exp3[exp3$gene_type %in% c("protein_coding"),]
lncRNA_exp2 <- exp3[exp3$gene_type %in% c("lncRNA"),]

# 最后去掉用于筛选的gene_type这一列：
mRNA_exp2 <- mRNA_exp2[,-2]
lncRNA_exp2 <- lncRNA_exp2[,-2]


# 两种方法结果一致
res <- identical(mRNA_exp,mRNA_exp2)
if(!isTRUE(res)){
  stop("mRNA 方法一和方法二结果不一致")
}
identical(lncRNA_exp,lncRNA_exp2)
if(!isTRUE(res)){
  stop("lncRNA 方法一和方法二结果不一致")
}

# 保存数据
data_path <- paste(project_path,"/data/",item,"/",item,".Rdata",sep = "")
save(
  # 做完合并的表达矩阵
  exp,
  # 添加gene_name列的表达矩阵
  exp2,
  # mRNA的表达矩阵
  mRNA_exp,
  #lncRNA的表达矩阵
  lncRNA_exp,
  file = c(data_path)
)