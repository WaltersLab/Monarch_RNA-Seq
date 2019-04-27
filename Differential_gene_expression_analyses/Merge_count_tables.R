#### Merge count tables procuded by STAR alinger for downstream DE analyses ####
# use "unstraned (column 2)" counts due to library type 

## Input STAR count files for all samples ----
file_list <- list.files(pattern = "_ReadsPerGene.out.tab")  # get all file names
name_list <- numeric(length(file_list))  # make a list for sample names 
list <- strsplit(file_list, split = "_")  # split 
for (i in 1:length(file_list)){  # extract names 
  name_list[i] <- list[[i]][1]
}
name_list <- gsub("-", "_", name_list)  # correct format 
# read in all dataframe as a list
all_data <- lapply(file_list, FUN = read.delim, header = F, as.is = T)  
sapply(all_data, dim)  # check dimensions

## Merge dataframes ----
# creat an empty dataframe 
DATA <- data.frame(matrix(NA, nrow = sapply(all_data, dim)[1] - 4, ncol = length(file_list)))
dim(DATA)   # nrow = all genes - headers, ncol = num samples
colnames(DATA) <- name_list  # add sample name 
# Fill in the dataframe 
for (i in 1:length(file_list)){
  data_in <- all_data[[i]][-(1:4), 1:2]
  if (i == 1){  # add row names 
    rownames(DATA) <- data_in$V1 
  }
  if (sum(data_in$V1 != rownames(DATA)) == 0){   # make sure gene names and orders all matched 
    DATA[, i] <- data_in$V2    
  }
  else{
    break
  }
}
head(DATA)

## Output 
write.table(DATA, file = "Merged_counts_all.txt")
