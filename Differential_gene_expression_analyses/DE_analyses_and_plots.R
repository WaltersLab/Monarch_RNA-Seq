#### Differential Gene Expression analysis ####
# DGE analyses were performed using Bioconductor package "EdgeR" 
# Generate MA plots, Volcano plots, and Heatmaps 
# Also test for specific gene sets of interest (immune and detoxification)
# Reference[1] https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# Reference[2] http://combine-australia.github.io/RNAseq-R/index.html
# Reference[3] http://www.nathalievilla.org/doc/html/solution-edgeR-rnaseq.html


### load libraries
library(edgeR)
library(RColorBrewer)
library(gplots)
library(limma)
library(grid)

### Import datasets----
# Input a merged expression count table generated previously
count_data <- read.table("Merged_counts_all.txt", header = T)
sample_id <- colnames(count_data)
# Input a list of sample information 
samples <- read.table("Sample_info_all.txt", header = T)  # input sample information 
samples <- samples[match(sample_id, samples$ID), ]  # sort the matrix according to the count dataframe
samples$ID == sample_id  # double checking 
colnames(count_data) <- samples$sample_name # replace IDs with sample name
# Create subsets
gut <- which(samples$tissue == "g") # gut tissue only 
body <- which(samples$tissue == "w")  # body tissue only 
# Input a list of monarch immune genes 
IG_list <- read.table("Immune_gene_list.txt", header = F)
# Input lists of monarch detoxification genes 
CYP <-  read.table("CYP_genelist.txt", header = F)   #Cytochrome P450 
UGT <-  read.table("UGT_genelist.txt", header = F)   # UDP-glycosyltransferase
ABC <-  read.table("ABC_genelist.txt", header = F)   # ABC transporter
GST <-  read.table("GST_genelist.txt", header = F)   # glutathione S-transferase

### DGE estimation and exploratory graphs ----
## All in a function 
DGE_estimation <- function(SUB, EXPLORE, outdir){  
# SUB = subset to use (gut or body);  EXPLORE = exploratory mode (T or F); outdir = output directory

## 1) Subset the data: either gut or body
data_sub <- count_data[, SUB]
samples_sub <- samples[SUB,] 
print("Sample info:")
print(samples_sub)

## 2) Create a data object (DGEList format) for edgeR analysis 
DGE <- DGEList(counts = data_sub)  # make a DGElist object 
# Create a grouping factor
groups <- as.factor(paste(samples_sub$treat, samples_sub$plant, sep = "_")) 

## 3) Filtering and normalizing data
print("Total gene counts per sample:")
print(apply(DGE$counts, 2, sum))
keep <- rowSums(cpm(DGE) > 0)  >= 2 # get index for genes with cpm > 0 in at least two samples
DGE <- DGE[keep,]  #filtering based on the above
# reset library sizes
DGE$samples$lib.size <- colSums(DGE$counts)
# Normalizing the data
DGE_norm <- calcNormFactors(DGE)   # calc. scaling factor by TMM method

## 4) Data Exploration
if (EXPLORE == T){  # if doing this option 
  setwd(paste0(current_dir, "/", outdir))
  ## Quality plots
    logcounts_un <- cpm(DGE,log = TRUE) # Get log2 CPM for unnormalized samples
    logcounts_nor <- cpm(DGE_norm,log = TRUE) # Get log2 CPM for normalized samples
  # Check distributions of unnormalized samples
  png(file="Quality_plot_before_normalization.png")
  boxplot(logcounts_un, xlab = "", ylab = expression(Log[2]("Counts per Million")), las = 2)
  abline(h = median(logcounts_un),col = "blue")  # median logCPM
  title("Boxplots of LogCPMs (unnormalized)")
  dev.off()
  # Check distributions of normalized samples
  png(file="Quality_plot_after_normalization.png")
  boxplot(logcounts_nor, xlab = "", ylab = expression(Log[2]("Counts per Million")), las = 2)
  abline(h = median(logcounts_nor),col = "blue")  # median logCPM
  title("Boxplots of LogCPMs (normalized)")
  dev.off()
  
  ## MDS plot
  colors <- c("red4", "firebrick1", "blue4", "dodgerblue1")
  png(file = "MDS_plot.png", width = 440, height = 480)
  plotMDS(DGE_norm, method="bcv", col = colors[as.factor(groups)],
          main = "Multidimensional scaling plot for samples",
          cex.lab = 1.3, cex.axis = 1.2, cex.main = 1.5, cex = 0.8)
  legend("topright", levels(groups), col = colors, pch=15, cex = 1.2)
  dev.off()
  
  ## Heatmap
  var_genes <- apply(logcounts_nor, 1, var)  # estimate var. for each row in the logcounts
  select_var <- names(sort(var_genes, decreasing=TRUE))[1:500] # Get the gene names for the top 500 most variable genes
  highly_variable_lcpm <- logcounts_nor[select_var,] # Subset the matrix
  color <- brewer.pal(11,"RdYlBu") 
  color_heat <- colorRampPalette(color)
  col.group <- colors[as.factor(groups)] # colors for groups
  # Plot the heatmap
  png(file = "heatmap.png", width = 500, height = 400)
  par(oma = c(0.5,3,0.5,0.5), xpd = T)
  heatmap.2(highly_variable_lcpm,col = rev(color_heat(50)),trace = "none", 
          main = "", margin = c(5,6), labRow = "",
          ColSideColors = col.group, scale = "row", 
          key.par = list(mgp = c(1.6, 0.5, 0), mar = c(3, 2.5, 3, 1), cex = 0.7, cex.lab = 1.3, cex.main = 1.2))
		  legend(0.2, 1.2, levels(groups), col = colors, pch=15, cex = 1.0, horiz = T)
  dev.off()
  par(oma = c(0, 0, 0, 0), xpd = F)
}


## 5) Estimating Dispersion using GLMs
# reate a design matrix
design.mat <- model.matrix(~ 0 + groups)
# estimate dispersion 
DGE_final <- estimateDisp(DGE_norm, design.mat)  # estimate common, trended, and tagwise dispersions
if (EXPLORE == T){  # if doing this option 
  png(file = "glm_dispersion.png")
  plotBCV(DGE_final, main = "Estimated Dispersion by GLM")
  dev.off()
}

## 6) Differential Expression
# fit glm with NB
fit <- glmFit(DGE_final, design.mat) # model fitting
# see comparisons
print("Comparisons:")
print(colnames(fit))

# Return the normalized DGE dataframe and the GLM-fit 
RESULTS <- list(groups = groups, DGE_norm = DGE_norm, DGE_final = DGE_final, fit = fit)
return(RESULTS)

# change back the wd
setwd(current_dir)
}


### DE analyses and graphs -----
## All in a function 
DE_analysis <- function(comp, outdir, FIGMODE, FIGLAB){
  # comp:  1 = inf in all, 2 = inf in inc, 3 = inf in cur, 4 = plants
  # outdir: output directory 
  # FIGMODE: 1 = MA+volcano plots; 2 = heatmap; F = N/A 
  # FIGLAB: figure labels
setwd(paste0(current_dir, "/", outdir))

## 1) Likihood-Ratio Tests
# comp:  1 = inf in all, 2 = inf in inc, 3 = inf in cur, 4 = plants
if (comp == 1){
  LRT <- glmLRT(fit, contrast = c(1,1,-1,-1))   # inf vs uninf (in both plants)
} else if (comp == 2){
  LRT <- glmLRT(fit, contrast = c(0,1,0,-1))   # inf vs uninf (in inc)
} else if (comp == 3){
  LRT <- glmLRT(fit, contrast = c(1,0,-1,0))  # inf vs uninf (in cur)
} else if (comp == 4){
  LRT <- glmLRT(fit, contrast = c(1,-1,1,-1))   # for cur vs inc (in both trts)
}  

## 2) Summary of num. DE genes 
DE <- decideTestsDGE(LRT, adjust.method = "BH", p.value = 0.05) 
print("Summary of differentially expressed genes:")
print(summary(DE))  # -1 = down-regulatedl 0 = non-diff; 1 = up-regulated

## 3) MA plot for DE genes with FDR < 0.05
detags <- rownames(DGE_final)[as.logical(DE)]  # DE genes
# customize the FC axis label 
if (comp == 4){
	FC_axis <- expression(paste("Log"[2], " (cur:inc)"))
	} else {
		FC_axis <- expression(paste("Log"[2], " (Inf:Uninf)"))
}
# plot 
if (FIGMODE == F){
  png(file = "smear_plot.png")
} 
if (FIGMODE != 2){
  par(mar = c(4,4.3,5,3))
  plotSmear(LRT, de.tags=detags,
            ylab = FC_axis,
            xlab = expression(paste("Log"[2], " average expression (CPM)")),
            cex.lab = 1.3, cex.axis = 1.2, 
            ylim = c(-25, 15)
            )  
  abline(h = c(-1, 1), col = "blue", lty = 5)  # line indicating +-1 fold change
  # Add title, grid line, figure number for graphic purposes 
  if (FIGMODE != F){
    mtext(substitute(bold(x), list(x = FIGLAB[1])), side = 3, adj = -0.25, line = 0.8, cex = 1.3)
	if (FIGLAB[1] %in% c("(A)", "(C)")){
	X <- grconvertX(17, "user", "ndc")
	grid.lines(x = X, y = c(0.01, 0.99), gp = gpar(col = "darkgray", lty = 2, lwd = 2))
  }
	if (FIGLAB[1] == "(A)"){
	title(expression(bold(underline("Infected vs. Uninfected"))), cex.main = 1.6, line = 3.5)
  } else if (FIGLAB[1] == "(B)"){
	title(expression(underline(paste(bolditalic("A. curassavica "), bold("vs. "), bolditalic("A. incarnata")))), cex.main = 1.6, line = 3.5)
	}
  }
}
if (FIGMODE == F){
  dev.off()
}

## 4) Volcano plot 
voc <- topTags(LRT, n = nrow(LRT))  # all genes
voc_color <- numeric(nrow(voc))  # colors 
for (i in 1:nrow(voc)){
  if(voc$table$logFC[i] >= 2 & voc$table$FDR[i] < 0.05){
    voc_color[i] = "red"
  } else if (voc$table$logFC[i] <= -2 & voc$table$FDR[i] < 0.05){
    voc_color[i] = "blue"
  } else {
    voc_color[i] = "black"
  }
}
# plot 
if (FIGMODE == F){
  png(file = "volcano_plot.png")
} 
if (FIGMODE != 2){
    par(mar = c(4,4.3,5,3))
    plot(voc$table$logFC, -log10(voc$table$FDR), pch=19, cex=0.3, col = voc_color,
         xlab = FC_axis,
         ylab = expression(paste("-Log"[10], " (P-value)"))
         , ylim = c(0, 5), xlim = c(-20, 15),
         cex.lab = 1.3, cex.axis = 1.2
    )
    abline(h = -log10(0.05), lty = 2)
    abline(v = -2, lty = 3)
    abline(v = 2, lty = 3)
	# Add figure number for graphic purposes 
    if (FIGMODE != F){
      mtext(substitute(bold(x), list(x = FIGLAB[2])), side = 3, adj = -0.25, line = 0.8, cex = 1.3)
  }
}
if (FIGMODE == F){
  dev.off()
}

## 5) HeatMap on only DE genes 
logcounts_nor <- cpm(DGE_norm,log=TRUE) # Get log2 CPM for normalized samples
var_genes <- apply(logcounts_nor, 1, var)  # estimate var. for each row in the logcounts
# Get the gene names for only the DE genes
if (length(detags) >= 250){
  select_var_DE <- intersect(names(sort(var_genes, decreasing = TRUE)), detags)[1:250]  # first 250 genes
} else {
  select_var_DE <- intersect(names(sort(var_genes, decreasing=TRUE)), detags)  # all DE genes
}
highly_variable_lcpm_DE <- logcounts_nor[select_var_DE,] # Subset the matrix
colors <- c("red4", "firebrick1", "blue4", "dodgerblue1")
color <- brewer.pal(11,"RdYlBu")  
color_heat <- colorRampPalette(color)
col.group <- colors[as.factor(groups)] # colors for groups
# Plot the heatmap
if (length(detags) > 10){   # no need to plot if having too few genes 
  if (FIGMODE == F){
    png(file="heatmap_DE.png", width = 500, height = 400)
  }
  if (FIGMODE != 1){
  par(oma = c(0.5,3,0.5,0.5), xpd = T)
  heatmap.2(highly_variable_lcpm_DE,col=rev(color_heat(50)),trace="none", 
          main="", margin = c(5,6), labRow = "",
          ColSideColors=col.group, scale="row", 
          key.par = list(mgp = c(1.6, 0.5, 0), mar = c(3, 2.5, 3, 1), cex = 0.7, cex.lab = 1.3, cex.main = 1.2))
	legend(0.3, 1.2, levels(groups), col = colors, pch=15, cex = 1.0, horiz = T)	  
  }
  # Add figure number for graphic purposes 
  if (FIGMODE == 2){
	  mtext(substitute(bold(x), list(x = FIGLAB[1])), side = 3, adj = -0.2, line = 3)
    }
  if (FIGMODE == F){
    dev.off()
  }
  par(oma = c(0, 0, 0, 0), xpd = F) 
}

## 6) Print out the top 15 up-reduated and down-regulated genes 
LRT_all <- topTags(LRT, n = nrow(count_data), sort.by = "PValue")  # all genes
sig_LRT_all <- as.data.frame(LRT_all[which(LRT_all$table$FDR < 0.05),]) # sig. DE genes
Pos_sig <- sig_LRT_all[which(sig_LRT_all$logFC > 0), ]  # up-reg
Neg_sig <- sig_LRT_all[which(sig_LRT_all$logFC < 0), ]  # down-reg
upreg_top_15 <- Pos_sig[order(Pos_sig$FDR, decreasing = F), ][1:15, ]  # top 15 up-reg
downreg_top_15 <- Neg_sig[order(Neg_sig$FDR, decreasing = F), ][1:15, ] # top 15 down-reg
print("Top 15 up-regulated genes:")
print(upreg_top_15)
print("Top 15 down-regulated genes:")
print(downreg_top_15)
#output
if (FIGMODE == F){
  write.table(rownames(Pos_sig), file = "DE_genes_upreg.txt",sep = "\t", row.names = F, col.names = F, quote = F)   # all sig. upreg. genes (for GO analysis)
  write.table(rownames(Neg_sig), file = "DE_genes_downreg.txt", sep = "\t", row.names = F, col.names = F, quote = F)   # all sig. downreg. genes (for GO analysis)
  write.table(upreg_top_15, file = "upreg_top_15.txt", sep = "\t")  # top 15 upreg.
  write.table(downreg_top_15, file = "downreg_top_15.txt", sep = "\t")  # top 15 downreg.
}

## 7) Find immune genes among all DE genes 
LRT_all <- topTags(LRT, n = nrow(count_data), sort.by = "PValue")  # all genes
sig_LRT_all <- as.data.frame(LRT_all[which(LRT_all$table$FDR < 0.05),]) # sig. DE genes
DE_genes <- rownames(sig_LRT_all)
DE_IG_genes <- sig_LRT_all[intersect(DE_genes, IG_list[,1]),]  # find intersections between the two 
print("DE immune genes:")
print(DE_IG_genes)
# output
if (FIGMODE == F){
  write.table(DE_IG_genes, file = "Immune_DE_list.txt", sep = "\t")
}

## 8) Find detoxification genes (only for plant comparisons)
if (comp == 4){
# expressed genes 
exp_genes <- rownames(DGE_norm$counts)  # defined as this 
print("Total expressed genes:")
print(length(exp_genes))  # no. of expressed genes

# Chose the 3 importnat canonical detox. genes Based on Bimbaum et al 2017 Mol Ecol 
# CYPs
print("Expressed CYP genes:")
print(length(intersect(exp_genes, CYP[,1]))) # expressed CYPs
sig_CYP <- sig_LRT_all[intersect(DE_genes, CYP[,1]),]  # find matchings 
print("Up-reg. CYP genes:")
print(length(which(sig_CYP$logFC > 0)))  # sig. up-reg. genes
print("Down-reg. CYP genes:")
print(length(which(sig_CYP$logFC < 0))) # sig. down-reg. genes

# UGTs
print("Expressed UGT genes:")
print(length(intersect(exp_genes, UGT[,1]))) # expressed UGTs
sig_UGT <- sig_LRT_all[intersect(DE_genes, UGT[,1]),]  # find matchings 
print("Up-reg. UGT genes:")
print(length(which(sig_UGT$logFC > 0)))  # sig. up-reg. genes
print("Down-reg. UGT genes:")
print(length(which(sig_UGT$logFC < 0))) # sig. down-reg. genes

# ABCs 
print("Expressed ABC genes:")
print(length(intersect(exp_genes, ABC[,1]))) # expressed ABCs
sig_ABC <- sig_LRT_all[intersect(DE_genes, ABC[,1]),]  # find matchings 
print("Up-reg. ABC genes:")
print(length(which(sig_ABC$logFC > 0)))  # sig. up-reg. genes
print("Down-reg. ABC genes:")
print(length(which(sig_ABC$logFC < 0))) # sig. down-reg. genes

# GSTs
print("Expressed GST genes:")
print(length(intersect(exp_genes, GST[,1]))) # expressed GSTs
sig_GST <- sig_LRT_all[intersect(DE_genes, GST[,1]),]  # find matchings 
print("Up-reg. GST genes:")
print(length(which(sig_GST$logFC > 0)))  # sig. up-reg. genes
print("Down-reg. GST genes:")
print(length(which(sig_GST$logFC < 0))) # sig. down-reg. genes

# outputs 
if (FIGMODE == F){
  write.table(sig_CYP, file = "CYP_DE_list.txt", sep = "\t")
  write.table(sig_UGT, file = "UGT_DE_list.txt", sep = "\t")
  write.table(sig_ABC, file = "ABC_DE_list.txt", sep = "\t")
  write.table(sig_GST, file = "GST_DE_list.txt", sep = "\t")
}
}

# change back the wd
setwd(current_dir)
}


### Run the function for gut samples ----
OUT <- DGE_estimation(SUB = gut, EXPLORE = T, outdir = "DE_results/Gut")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
DE_analysis(comp = 1, outdir = "DE_results/Gut/Inf_in_all", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 2, outdir = "DE_results/Gut/Inf_in_INC", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 3, outdir = "DE_results/Gut/Inf_in_CUR", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 4, outdir = "DE_results/Gut/Plant_in_all", FIGMODE = F, FIGLAB = "")

### Run the function for body samples ----
OUT <- DGE_estimation(SUB = body, EXPLORE = T, outdir = "DE_results/Body")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
DE_analysis(comp = 1, outdir = "DE_results/Body/Inf_in_all", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 2, outdir = "DE_results/Body/Inf_in_INC", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 3, outdir = "DE_results/Body/Inf_in_CUR", FIGMODE = F, FIGLAB = "")
DE_analysis(comp = 4, outdir = "DE_results/Body/Plant_in_all", FIGMODE = F, FIGLAB = "")

### Generate Figures ---- 
## Fig. 3 (MA plots and volcano plots for gut samples)
tiff(file = "FIGURE3.tif", width = 2400, height = 2400, res = 300)
OUT <- DGE_estimation(SUB = gut, EXPLORE = F, outdir = "")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
layout(matrix(c(1,2,3,4), 2, 2))
DE_analysis(comp = 1, outdir = "", FIGMODE = 1, FIGLAB = c("(A)", "(C)"))
DE_analysis(comp = 4, outdir = "", FIGMODE = 1, FIGLAB = c("(B)", "(D)"))
dev.off()

## Fig. 4 (MA plots and volcano plots for body samples)
tiff(file = "FIGURE4.tif", width = 2400, height = 2400, res = 300)
OUT <- DGE_estimation(SUB = body, EXPLORE = F, outdir = "")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
layout(matrix(c(1,2,3,4), 2, 2))
DE_analysis(comp = 1, outdir = "", FIGMODE = 1, FIGLAB = c("(A)", "(C)"))
DE_analysis(comp = 4, outdir = "", FIGMODE = 1, FIGLAB = c("(B)", "(D)"))
dev.off()

## Fig. 5 (heatmaps for gut and body samples)
tiff(file = "FIGURE5A.tif", width = 2400, height = 1600, res = 300)
OUT <- DGE_estimation(SUB = gut, EXPLORE = F, outdir = "")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
DE_analysis(comp = 4, outdir = "", FIGMODE = 2, FIGLAB = "(A)")
dev.off()

tiff(file = "FIGURE5B.tif", width = 2400, height = 1600, res = 300)
OUT <- DGE_estimation(SUB = body, EXPLORE = F, outdir = "")
groups <- OUT$groups
DGE_norm <- OUT$DGE_norm
DGE_final <- OUT$DGE_final
fit <- OUT$fit
DE_analysis(comp = 4, outdir = "", FIGMODE = 2, FIGLAB = "(B)")
dev.off()

