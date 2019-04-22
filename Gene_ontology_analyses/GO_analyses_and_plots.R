#### Gene Ontology (GO) Enrichcment Analyses ####
# Monarch GO annotations previously generated using "Pannzer2" (http://ekhidna2.biocenter.helsinki.fi/sanspanz/). 
# A custom organism package built ioconductor package using "AnnotationForge".
# Enrichment tests were performed using Bioconductor package "clusterProfiler".
# Reference[1] https://bioconductor.org/packages/release/bioc/vignettes/AnnotationForge/inst/doc/MakingNewOrganismPackages.html
# Reference[2] https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

### load libraries
library(AnnotationForge)
library(clusterProfiler)
library(ggplot2)
library(gridExtra)
library(grid)


### Step 1: Build an "OrgDb" object for monarchs using AnnotationForge---- 
## Import datasets 
# GO annotations 
Dp_GO_raw <- read.delim("GO.out.txt", header = T)  # input the annotation output from PANNZER2
str(Dp_GO_raw)
# modify the dataframe structure and formats 
Dp_GO <- cbind(Dp_GO_raw[, c(1,3)], "IEA")  # extract things needed and add a column of evidance code (all "IEA" here)
Dp_GO <- Dp_GO[Dp_GO[, 2] != "", ]  # remove NA entries (if any)
colnames(Dp_GO) <- c("GID", "GO", "EVIDENCE")  # add headers
Dp_GO$GID <- gsub("-PA", "", Dp_GO$GID)  # modify gene ID format (no "-PA")
Dp_GO$GO <- paste0("GO:", sprintf("%07d", Dp_GO$GO)) # modify GO ID format (add GO: and zeros)
Dp_GO$GO <- as.factor(Dp_GO$GO)  #convert type 
Dp_GO$GID <- as.character(Dp_GO$GID)
str(Dp_GO)

# All monarch genes
Dp_genes_raw <- read.table("Dp_gene_info.txt", header = T, sep = "\t") # input a dataframe of all monarch genes and their location information (modified from data of https://doi.org/10.1534/g3.117.300187)
str(Dp_genes_raw)
Dp_genes <- Dp_genes_raw[, c(2, 1)]  # extract 
Dp_genes <- Dp_genes[which(Dp_genes[, 2] != ""), ]  # remove NA entries (if any)
colnames(Dp_genes) <- c("GID", "CHROMOSOME")  # add header. note that the chromosomes are just scaffolds. 
Dp_genes$GID <- as.character(Dp_genes$GID)   # convert type 
Dp_genes <- Dp_genes[order(Dp_genes$GID), ] # sort 
str(Dp_genes)

# Check coverage 
length(unique(Dp_GO$GID))  # genes that have GO annotations 
length(unique(Dp_GO$GID))/nrow(Dp_genes)  # ratio 

# Use the "makeOrgPackage" function to create a custom package
makeOrgPackage(chromosome=Dp_genes, go=Dp_GO,
               version="0.1",
               maintainer="WH Tan <wtan4@emory.edu>",
               author="WH Tan <wtan4@emory.edu>",
               outputDir = "C:/R/Dp_GO",
               tax_id="278856",
               genus="Danaus",
               species="plexippus",
               goTable="go")

# Install the custom package
install.packages("C:/R/Dp_GO/org.Dplexippus.eg.db", repos=NULL, type = "source")  # type=source is needed for Windows


### Step 2: Import the RNA-seq results ---- 
# Import lists of significantly differentially expressed (DE) genes obtained from the DE analyses
DE_list <- list.files(pattern = "DE_genes_")  # list all DE result files 
all_data <- lapply(DE_list, read.table)  # read them in 
for (i in 1:length(DE_list)){  # assign each of them and convert name and type 
  assign(gsub("DE_genes_|.txt", "", DE_list[i]), as.character(all_data[[i]]$V1))
}

# Import a dataframe of all expressed genes 
RNASEQ_table <- read.table("Merged_counts_all.txt")   # all raw data 
str(RNASEQ_table)
RNASEQ_univ <- rownames(RNASEQ_table[which(rowSums(RNASEQ_table) > 0), ]) # extract the gene IDs of genes that have expressed in at least one individual 
str(RNASEQ_univ)


### Step 3: Load the monarch package for enrichment tests ----  
library(org.Dplexippus.eg.db)  # the monarch package 


### Step 4: Enrichment analyses ---- 
## Make a function for the enrichment test using ClusterProfiler
enrichment_test <- function(subset){   # subset  = the group of interest 
# Set the group of interest and check numbers 
gene_list <- subset  # set the group of interest 
N <- length(gene_list)  # see how many genes 
M <- length(intersect(gene_list, unique(Dp_GO$GID))) # see how many have GO annotations 
print("Num. of genes, Num. of genes with GO annotations, the ratio:")
print(c(N, M, M/N)) # see N, M, and the ratio 

# Run the over-representation test
ego <- enrichGO(gene = gene_list,     # group of (DE) genes of interest 
                universe = as.character(RNASEQ_univ),   # the background set (all genes expressed in the RNAseq)
                OrgDb = org.Dplexippus.eg.db,   
                keyType = "GID", 
                ont = "ALL", 
                pAdjustMethod = "BH",  
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, 
                readable = FALSE)
return(ego)
}

## Run the function for all groups of interest
ego_upreg_gut <- enrichment_test(upreg_gut)
ego_downreg_gut <- enrichment_test(downreg_gut)
ego_upreg_body <- enrichment_test(upreg_body)
ego_downreg_body <- enrichment_test(downreg_body)

## output test results 
out1 <- head(ego_upreg_gut, 14849)
write.table(out1, file = "upreg_in_gut.txt", sep = "\t")
out2 <- head(ego_downreg_gut, 14849)
write.table(out2, file = "downreg_in_gut.txt", sep = "\t")
out3 <- head(ego_upreg_body, 14849)
write.table(out3, file = "upreg_in_body.txt", sep = "\t")
out4 <- head(ego_downreg_body, 14849)
write.table(out4, file = "downreg_in_body.txt", sep = "\t")


### Step 5: Plotting ----
P1 <- dotplot(ego_upreg_gut, showCategory = 20, font.size = 12)
P2 <- dotplot(ego_downreg_gut, showCategory = 20, font.size = 12)
P3 <- dotplot(ego_upreg_body, showCategory = 20, font.size = 12)
P4 <- dotplot(ego_downreg_body, showCategory = 20, font.size = 12)

# Fig. 6: the gut subset 
tiff("FIGURE6.tif", width = 2400, height = 3600, res = 300)
grid.arrange(P1, P2, nrow = 2) 
grid.text("(A)", x = unit(0.02, "npc"), y = unit(0.98, "npc"), gp=gpar(fontsize = 17, fontface = "bold")) 
grid.text("(B)", x = unit(0.02, "npc"), y = unit(0.48, "npc"), gp=gpar(fontsize = 17, fontface = "bold")) 
dev.off()

# Fig. 7: the body subset 
tiff("FIGURE7.tif", width = 2400, height = 2400, res = 300)
grid.arrange(P3, P4, heights=c(0.3, 0.7), nrow = 2)
grid.text("(A)", x = unit(0.02, "npc"), y = unit(0.98, "npc"), gp=gpar(fontsize = 17, fontface = "bold"))
grid.text("(B)", x = unit(0.02, "npc"), y = unit(0.48, "npc"), gp=gpar(fontsize = 17, fontface = "bold"))
dev.off()
