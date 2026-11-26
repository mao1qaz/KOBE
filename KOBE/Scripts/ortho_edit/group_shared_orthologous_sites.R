#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2){
  stop("Usage: Rscript info_group_shared.R input.csv output.csv")
}

input_file <- args[1]
output_file <- args[2]

#-------------------------------
# 读取数据
#-------------------------------
df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

#-------------------------------
# 自动识别物种
#-------------------------------
cds_cols <- grep("_cds_pos$", colnames(df), value = TRUE)
species <- sub("_cds_pos$", "", cds_cols)

cat("Detected species:\n")
print(species)
cat("Number of species:", length(species), "\n")

#-------------------------------
# 生成 ID2 列
#-------------------------------
df$ID2 <- "unfound"

for(i in 1:nrow(df)){
  id_prefix <- sub(":.*", "", df$ID[i])
  
  for(col in cds_cols){
    val <- df[i, col]
    if(is.na(val)) next
    gene_name <- sub(":.*", "", val)
    if(gene_name == id_prefix){
      df$ID2[i] <- val
      break
    }
  }
}

#-------------------------------
# 检查 ID2 是否有重复
#-------------------------------
id2_valid <- df$ID2[df$ID2 != "unfound"]

if(any(duplicated(id2_valid))){
  cat("There are duplicated ID2 values.\n")
} else {
  cat("No duplicated ID2 values found.\n")
}

#-------------------------------
# 生成 Shared 列（T/F）
#-------------------------------
df$Shared <- "F"
all_ID2 <- df$ID2

for(i in 1:nrow(df)){
  cds_values <- unlist(df[i, cds_cols])
  cds_values <- cds_values[cds_values != "NOT_FOUND:0"]
  if(any(cds_values %in% all_ID2[-i])){
    df$Shared[i] <- "T"
  }
}

#-------------------------------
# 生成 ID2_group 列（排序+去重）
#-------------------------------
df$ID2_group <- vector("list", nrow(df))

for(i in 1:nrow(df)){
  cds_values <- unlist(df[i, cds_cols])
  cds_values <- cds_values[cds_values != "NOT_FOUND:0"]
  shared_id2 <- cds_values[cds_values %in% all_ID2[-i]]
  df$ID2_group[[i]] <- sort(unique(c(df$ID2[i], shared_id2)))
}

#-------------------------------
# 去重规则（保留 ID2_group 多的）
#-------------------------------
keep <- rep(TRUE, nrow(df))

for(i in 1:(nrow(df)-1)){
  if(!keep[i]) next
  for(j in (i+1):nrow(df)){
    if(!keep[j]) next
    set_i <- df$ID2_group[[i]]
    set_j <- df$ID2_group[[j]]
    
    # 完全相同
    if(setequal(set_i, set_j)){
      keep[j] <- FALSE
    } else if(all(set_i %in% set_j)){ # i被包含
      keep[i] <- FALSE
      break
    } else if(all(set_j %in% set_i)){ # j被包含
      keep[j] <- FALSE
    }
  }
}

df_final <- df[keep, ]

#-------------------------------
# 输出 CSV
#-------------------------------
df_final$ID2_group <- sapply(df_final$ID2_group, function(x) paste(x, collapse=";"))

write.csv(df_final, output_file, row.names = FALSE, quote = FALSE)
cat("Finished. Output written to:", output_file, "\n")