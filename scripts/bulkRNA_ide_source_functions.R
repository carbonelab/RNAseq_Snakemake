

# Functions for data exploration of bulk RNAseq data
# 1. clean_gene_info: read in a gene info file and produce clean colnames, remove dups
# 2. compile_readcounts: create counts table from directory of STAR ReadsPerGene.tab files
# 3. plot_total_counts: create barplot of total raw counts per sample
# 4. plot_top_count_genes: create barplot of top n gene by average count across all samples
# 5. lcf_edger: low count filter via edgeR method
# 6. var_stable_transform: perform variance stabilizing transformation of counts for visualization
# 7. sample_dist_heatmap: plot a Euclidean distance sample correlation heatmap
# 8. deseq2_pca: plot a simple PCA biplot of PC1 vs PC2 via DESeq2's plotPCA function

#----------------------------------------------------------------------------------------------------#

#############
# LIBRARIES #
#############

library(parallel)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(stringr)
library(readxl)
library(ggrepel)

#----------------------------------------------------------------------------------------------------#

############----------------------------------------------------------------------------------#
# FUNCTION # clean_gene_info: to read in a gene info txt file and clean up for downstream use #
############----------------------------------------------------------------------------------#
# INPUTS: gene_info: path to gene info txt file (tab-separated)
#         gene_id_col: name of column in file representing GeneID
#         gene_name_col: name of column in file representing gene name
#         description: whether or not a gene description column is present
#         gene_desc: name of column in file representing gene description
#         new_id: desired name of GeneID column
#         new_name: desired name of gene name column
#         new_desc: desired name of gene description column
# OUTPUTS: a cleaned gene_info data.frame
clean_gene_info <- function(gene_info, gene_id_col="GeneID", gene_name_col="gene.name", description=NULL, gene_desc="gene.description", new_id="GeneID", new_name="gene.name", new_desc="gene.description") {
  # read in gene info file and store as df
  ref <- read.delim(gene_info, header=TRUE)

  # remove any duplicate GeneID rows
  ref <- ref[!duplicated(ref[[gene_id_col]]), ]
  # change GeneID colname
  colnames(ref)[colnames(ref)==gene_id_col] <- new_id
  # change gene.name colname to "gene.name"
  colnames(ref)[colnames(ref)==gene_name_col] <- new_name

  if (!(is.null(description))) {
    # change gene description colname
    colnames(ref)[colnames(ref)==gene_desc] <- new_desc
  }

  return(ref)
}

#-----------------------------------------------------------------------------------------------#

############--------------------------------------------#
# FUNCTION # compile_readcounts: to create counts table #
############--------------------------------------------#
# INPUTS: dir: path to directory containing STAR ReadsPerGene.out.tab files
#         strand: number signifying which column of counts to keep
#         metadat: path to metadata excel file, first column is Sample_Name
#         outfile: name of output counts table text file
#         suffix: ending of ReadsPerGene files to grab (and to remove)
# OUTPUTS: 1. a raw counts table tab-separated text file
#          2. returns a raw counts data.frame object
compile_readcounts <- function(dir, strand, metadat, outfile="raw_counts.txt", suffix=".star.ReadsPerGene.out.tab") {
  # create a list that contains the path to all input files
  counts_list <- list.files(dir, paste0(suffix, "$"), full=T)
  names(counts_list) <- counts_list
  names(counts_list) <- gsub(suffix,'',names(counts_list))
  names(counts_list) <- gsub(dir, '',names(counts_list))
  names(counts_list) <- gsub("/", '',names(counts_list))

  # replace the current names (the fastq prefix) with desired sample names from metadata excel file
  # read in the metadata excel file
  my_meta <- as.data.frame(read_excel(metadat))
  colnames(my_meta)[1] <- "sample"
  # create vector of sample names from matching fastq with sample in metadata file
  my_samples <- as.character(sapply(names(counts_list), function(x) my_meta[my_meta$fastq == x, "sample"]))
  names(counts_list) <- my_samples

  # read in each counts file as a data frame and store in a list of df's
  counts_dfs <- mclapply(counts_list, function(x) read.table(x), mc.cores=8)

  # convert the correct counts column to numeric
  for (i in 1:length(counts_dfs)) {counts_dfs[[i]][,strand] <- as.numeric(counts_dfs[[i]][,strand])}

  # change the name of the correct column to be the sample name
  for (i in 1:length(counts_dfs)) {colnames(counts_dfs[[i]])[strand] <- names(counts_dfs)[i]}

  # change the name of column 1 to "GeneID"
  for (i in 1:length(counts_dfs)) {colnames(counts_dfs[[i]])[1] <- "GeneID"}

  # drop unnecessary counts columns
  for (i in 1:length(counts_dfs)) {counts_dfs[[i]] <- subset(counts_dfs[[i]], select = colnames(counts_dfs[[i]][c(1,strand)]))}

  # remove rows 1,2,3,4
  for (i in 1:length(counts_dfs)) {counts_dfs[[i]] <- counts_dfs[[i]][-c(1,2,3,4), ]}

  # merge all data frames on "GeneID"
  my_table <- Reduce(function(df1, df2) merge(df1, df2, by = "GeneID", all = TRUE), counts_dfs)

  # The order of the sample columns in the counts table is currently in the order of the
  #  ReadPerGene files from the star directory
  # rearrange sample columns of raw counts table to be in same order as rows of metadata table
  # if the desired sample names in metadata table have dashes, they will have been replaced in colnames of raw counts table as periods
  # replace the dashes with periods in the first column of metadata table
  #my_meta$sample <- gsub("-", ".", my_meta$sample)
  rownames(my_table) <- my_table[ ,1]
  my_table[ ,1] <- NULL
  my_table2 <- my_table[ ,match(my_meta$sample, colnames(my_table))]
  my_table2$GeneID <- rownames(my_table2)
  my_table2 <- my_table2[ ,c(ncol(my_table2),1:ncol(my_table2)-1)]

  # export the final table
  write.table(my_table2, outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  return(my_table2)
}

#----------------------------------------------------------------------------------------------#

############------------------------------------------------------------------#
# FUNCTION # plot_total_counts: create barplot of total raw counts per sample #
############------------------------------------------------------------------#
# INPUTS: raw_counts: raw counts data.frame (only samples as columns, GeneID is rownames)
#         outfile_prefix: desired name of plot file (excluding the suffix like .png or .pdf)
#         plot_type: type of plotfile, ex: "png" or "pdf"
# OUTPUTS: 1. a barplot of total counts file
#          2. returns the myplot ggplot2 object
plot_total_counts <- function(raw_counts=raw_counts, outfile_prefix="total_counts_barplot", plot_type="png") {
  # make a df of total counts per sample and then plot (refer to coldata for variables)
  tcounts.df <- data.frame(sample=colnames(raw_counts))
  tcounts.df$total_counts <- colSums(raw_counts)

  # create barplot of total counts
  myplot <- ggplot(tcounts.df, aes(x=sample,y=total_counts, fill=sample)) +
    labs(x="Sample", y="Total Gene Counts", title="Total Gene Counts per Sample after Alignment") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_bar(stat="identity") +
    theme_minimal()+
    theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #ggsave(filename=paste("total_counts_barplot", plot_type, sep="."), plot=myplot)

  return(myplot)
}

#----------------------------------------------------------------------------------------------------------#

############--------------------------------------------------------------------#
# FUNCTION # plot_top_count_genes: to plot the top genes based on average count #
############--------------------------------------------------------------------#
# PURPOSE: produce a barplot of the top genes by average count across all samples
#          export a barplot file and a corresponding text file
#          if gene info is provided, x axis labels will be gene names
#            - if a GeneID has no gene name, the label will be the GeneID
#            - if no gene info is provided, x axis labels will be GeneID
# INPUTS: raw_counts: raw counts data.frame (only samples as columns, GeneID is rownames)
#         n: number of top genes to plot, default 20
#         gene_info: (optional) gene info data frame (output of clean_gene_info function, must have GeneID and gene.name cols)
#         outfile: name of output text file (n genes ordered by avg count)
#         plotfile_prefix: prefix of plot file
#         plot_type: type of plot file ("png" or "pdf")
# OUTPUTS: 1. barplot file of top genes
#          2. max_counts.txt file, n genes ordered by avg count
#          2. returns ggplot2 object
plot_top_count_genes <- function(raw_counts=raw_counts, n=20, gene_info=NULL, outfile="max_counts.txt", plotfile_prefix="max_counts_barplot", plot_type="png") {
  # Add an average column to get the average counts for each gene
  raw_counts$avg <- rowMeans(raw_counts)
  # Re-add GeneID column
  raw_counts$GeneID <- rownames(raw_counts)
  # re-order to place GeneID first
  raw_counts <- raw_counts[ ,c(ncol(raw_counts),1:ncol(raw_counts)-1)]

  # get the top n most abundant genes
  max_counts <- raw_counts[order(raw_counts$avg, decreasing=T)[1:n], ]

  # if gene_info is provided, use it to get gene names
  if (!(is.null(gene_info))) {
    # merge max_counts table with subsetted ref table on GeneID
    max_counts_final <- merge(max_counts, gene_info, by="GeneID")

    # reorder rows by avg count
    max_counts_final <- max_counts_final[order(max_counts_final$avg, decreasing=T), ]

    # fill missing gene.name values with GeneID
    max_counts_final$gene.name <- ifelse(max_counts_final$gene.name == "", max_counts_final$GeneID, as.character(max_counts_final$gene.name))
  
    # specify the order of the gene.name factor levels for correct x-axis order on plot
    max_counts_final$gene.name <- factor(max_counts_final$gene.name, levels=unique(max_counts_final$gene.name))

    # create barplot of top Genes with Max Counts on Avg
    myplot <- ggplot(max_counts_final, aes(x=gene.name,y=avg)) +
      labs(x="Gene", y="Avg Gene Count", title="Genes with Highest Average Raw Counts Across All Samples") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(legend.position="none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #ggsave(filename=paste(plotfile_prefix, plot_type, sep="."), plot=myplot)
  }

  else if (is.null(gene_info)) {
    max_counts_final <- max_counts

    # reorder rows by avg count
    max_counts_final <- max_counts_final[order(max_counts_final$avg, decreasing=T), ]
    
    # specify the order of the gene.name factor levels for correct x-axis order on plot
    max_counts_final$GeneID <- factor(max_counts_final$GeneID, levels=unique(max_counts_final$GeneID))

    # create barplot of top Genes with Max Counts on Avg
    myplot <- ggplot(max_counts_final, aes(x=GeneID, y=avg)) +
      labs(x="Gene", y="Avg Gene Count", title="Genes with Highest Average Raw Counts Across All Samples") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      geom_bar(stat="identity") +
      theme_minimal() +
      theme(legend.position="none") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    #ggsave(filename=paste(plotfile_prefix, plot_type, sep="."), plot=myplot)
  }

    # re-order columns
    #max_counts_final <- max_counts_final[ ,c(1,(ncol(max_counts_final) - 1), ncol(max_counts_final), (ncol(max_counts_final)-2), 2:(ncol(max_counts_final) -3))]
    # reorder rows by avg count
    #max_counts_final <- max_counts_final[order(max_counts_final$avg, decreasing=T), ]

    # export max counts table
    #write.table(max_counts_final, outfile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  # specify the order of the gene.name factor levels for correct x-axis order on plot
  #max_counts_final$gene.name <- factor(max_counts_final$gene.name, levels=unique(max_counts_final$gene.name))

  # create barplot of top Genes with Max Counts on Avg
  #myplot <- ggplot(max_counts_final, aes(x=gene.name,y=avg)) +
    #labs(x="Gene Name", y="Avg Gene Count", title="Genes with Highest Average Raw Counts Across All Samples") +
    #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    #geom_bar(stat="identity") +
    #theme_minimal() +
    #theme(legend.position="none") +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #ggsave(filename=paste(plotfile_prefix, plot_type, sep="."), plot=myplot)

  return(myplot)
}

#-----------------------------------------------------------------------------------------------------------------#

############------------------------------------------#
# FUNCTION # lcf_edger: to low count filter via edgeR #
############------------------------------------------#
# INPUTS: raw_counts: raw_counts data.frame (samples only as columns)
#         group: vector denoting sample group membership
# OUTPUTS: 1. filtered_counts data.frame
lcf_edger <- function(raw_counts, group) {

  # create the DGEList data class from edgeR using the raw counts and group vector
  y <- DGEList(counts=raw_counts, group=group)

  # Filter low count genes
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]

  filtered_counts <- as.data.frame(y$counts)

  print(paste0("Removed ", nrow(raw_counts) - nrow(filtered_counts), " genes with low counts across samples"))
  print(paste0("Number of genes remaining: ", nrow(filtered_counts)))
  
  return(filtered_counts)
}

#-----------------------------------------------------------------------------------------------------------------#

############--------------------------------------------------------------------------------#
# FUNCTION # var_stable_transform: obtain transformed count data suitable for visualization #
############--------------------------------------------------------------------------------#
# PURPOSE: transform RNAseq count data for visualizations (remove dependence of variance on the mean)
#          Multiple methods exist. DESeq2 provides VST and rlog functions. edgeR uses logCPM.
#          Invoking a DESeq2 method will create the dds object, normalize for seq depth via estimateSizeFactors, and transform.
#          Invoking the edgeR method will create the DGEList object, normalize for seq depth via calcNormFactors, and transform.
# INPUTS: counts: a counts data.frame (usually filtered_counts), samples as columns
#         coldata: a coldata data.frame made from the metadata
#         method: one of "vst", "rlog", "logcpm"
#         formula_vars: desired variables for design creation (in order)
#                       only required if method is "vst" or "rlog"
#                       ex: formula_vars=c(0, "Group")
#         blind: whether or not blind is TRUE or FALSE for vst or rlog transform
#                only required if method is "vst" or "rlog"
#                blind set to FALSE is generally safe
#         group: group of interest for DGEList object creation
#                must be colname from coldata (often it's "Group")
#                only required if method is "logcpm"
# OUTPUTS: depends on the transformation method chosen. Either:
#          vst: a DESeqTransform object (genes as rows by samples as columns)
#          rlog: a DESeqTransform object (genes as rows by samples as columns)
#          logcpm: a data.frame of logcpm values
var_stable_transform <- function(counts, coldata=coldata, method, formula_vars=NULL, blind=NULL, group=NULL) {
  if (method %in% c("vst", "rlog")) {
    # Create the design from the formula variables
    design <- as.formula(paste("", paste(formula_vars, collapse= " + "), sep=" ~ "))
    # Create a dds object via DESeq2
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=design)
    dds <- estimateSizeFactors(dds)
    # Transform to stabilize variance
    if (method=="vst") {
      print("Performing vst transformation via DESeq2")
      print(paste0("Design formula specified as: ~ ", design[2]))
      vsd <- vst(dds, blind=blind)
      return(vsd)
    }
    else if (method=="rlog") {
      print("Performing rlog transformation via DESeq2")
      print(paste0("Design formula specified as: ~ ", design[2]))
      rlog <- rlog(dds, blind=blind)
      return(rlog)
    }
  }
  else if (method=="logcpm") {
    print("Performing logCPM transformation via edgeR")
    print(paste0("Group of interest specified as: ", group))
    # Create DGEList object via edgeR
    y <- DGEList(counts=counts, group=coldata[[group]])
    # Normalize (TMM normalization for RNA composition)
    y <- calcNormFactors(y, method = "TMM")
    # Log CPM Transform
    cpm.log <- as.data.frame(cpm(y, log=TRUE))
    return(cpm.log)
  }
}

#-------------------------------------------------------------------------------------------------------------------------#

############-----------------------------------------------------------------#
# FUNCTION # sample_dist_heatmap: plot a sample distance correlation heatmap #
############-----------------------------------------------------------------#
# INPUTS: norm_counts: normalized counts, output of var_stable_transform function
#                      note that it needs to be properly formatted (ex: t(assay(vsd)))
#         outfile: prefix of output plot file, default: eucl_dist_heatmap
#         plot_type: type of plotfile ("png" or "pdf")
#         print_to_screen: whether or not to print to the screen (such as in rmarkdown doc)
# OUTPUTS: a png or pdf heatmap plot file
#          optionally, prints pheatmap plot to screen
sample_dist_heatmap <- function(norm_counts, outfile="eucl_dist_heatmap", plot_type=plot_type, print_to_screen=NULL, plot_prefix="./") {
  # use "dist" function to calculate Euclidian distance between samples
  sampleDists <- dist(norm_counts)

  # visualize the distances in a heatmap
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

  if (plot_type=="png") {
    png(paste0(plot_prefix, outfile, ".", plot_type))
    pheatmap(sampleDistMatrix,
         main="Sample Distance Heatmap",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
    dev.off()
  } else if (plot_type=="pdf") {
    pdf(paste0(plot_prefix, outfile, ".", plot_type))
    pheatmap(sampleDistMatrix,
         main="Sample Distance Heatmap",
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
    dev.off()
  }

  if (!(is.null(print_to_screen))) {
    pheatmap(sampleDistMatrix,
             main="Sample Distance Heatmap",
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
  }	   
}

#----------------------------------------------------------------------------------------------------------------------------#

############---------------------------------------------#
# FUNCTION # deseq2_pca: use DESeq2 to plot a PCA biplot # 
############---------------------------------------------#
# INPUTS: norm_counts: normalized counts, must be either vst or rlog output of var_stable_transform function
#         intgroup: variable (colname of coldata) to be used for color of PCA points            
#                   default, "Group"
#         plot_type: type of plotfile ("png" or "pdf")
#         outfile: prefix of output plot file, default: deseq2_pca
# OUTPUTS: 
deseq2_pca <- function(norm_counts, intgroup="Group", outfile="deseq2_pca", plot_type=plot_type) {
  # use all genes passing low count filter for PCA
  dataVST <- plotPCA(norm_counts, intgroup=intgroup, returnData=TRUE)
  percentVar <- round(100 * attr(dataVST, "percentVar"))
  
  # Make PCA Plot
  myplot <- ggplot(dataVST, aes_string("PC1", "PC2", color=intgroup)) +
              geom_point(size=3) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance")) +
              ggtitle("Principal Component Analysis") +
              geom_text(aes(label=paste(name, sep="")),hjust=1, vjust=-.5) +
              theme_minimal()
  #ggsave(filename=paste(outfile, plot_type, sep="."), plot=myplot)

  return(myplot)
}

#---------------------------------------------------------------------------------------------------------------------------#

############----------------------------------#
# FUNCTION # function to describe pca results #
############----------------------------------#
# INPUTS: mypca: result of "prcomp" function
#
# OUTPUTS: sentences describing PCA results (number of PCs, % variance of PC1, etc)
#          returns num_pc
describe_pca <- function(mypca) {
  # Determine how many PC's were returned
  print(paste0("Number of PC's returned: ", ncol(mypca$x)))

  # Obtain the eigenvalues (can get values proportional to eigenvalues by taking sd^2)
  eigs <- mypca$sdev^2

  # Determine number of PC's with eigenvalue > 1 (considered important)
  print(paste0("Number of PC's with eigenvalue > 1: ", sum(eigs > 1)))

  # Determine how much of total variance is explained by first PC
  print(paste0("Percent of total variance explained by first PC: ", round((eigs[1]/sum(eigs))*100, digits=2), "%"))

  # How many PC's are needed to explain at least 80% of total variance
  my_sum = 0
  num_pc = 0
  for (i in 1:ncol(mypca$x)) {
    my_sum = my_sum + (eigs[i] / sum(eigs))
    num_pc = num_pc + 1
    if (my_sum >= 0.8) {
      break
    }
  }
  print(paste0("Number of PC's required to explain at least 80% of the variance: ", num_pc))

  return(num_pc)
}

#---------------------------------------------------------------------------------------------------------------------------#

############---------------------------------#
# FUNCTION # manually perform pca via prcomp #
############---------------------------------#
# INPUTS:  vsd: a variance stabilized object (vsd, rlog, etc)
#          n: number of features to use (default 500)
# OUTPUTS: produces a screeplot
#          describes the pca via custom function
#          returns the prcomp object mypca
perform_manual_pca <- function(vsd, n=500) {
  # select the features to use (highest variance)
  ntop = n
  Pvars <- rowVars(assay(vsd))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]

  # perform pca
  mypca <- prcomp(t(assay(vsd)[select,]))

  # describe pca
  describe_pca(mypca)

  # produce screeplot
  png("pca_screeplot.png")
  screeplot(mypca, npcs=length(mypca$sdev))
  dev.off()

  # get relative contribututions of each gene for each PC
  #aload <- abs(mypca$rotation) ## save absolute values
  #contr <- sweep(aload, 2, colSums(aload), "/")

  return(mypca)
}

#---------------------------------------------------------------------------------------------------------------------------#

############----------------------------------#
# FUNCTION # general custom plot pca function #
############----------------------------------#
# INPUTS: mypca: result of prcomp function
#         coldata: df of columns/metadata to use in PCA plot
#         pc_a: which PC to use (ex: "PC1")
#         pc_b: which PC to use (ex: "PC2")
#         color_var: column in coldata to color points by
#         label_var: column in coldata to label points
# 
# OUTPUTS: returns the ggplot object to be saved outside the function
plot_pca <- function(mypca, coldata, pc_a="PC1", pc_b="PC2", color_var, label_var="Sample_Name") {
  # grab the principle components as a new df
  my_pca_df <- as.data.frame(mypca$x)

  # add the relevant variable columns from coldata to the pca df
  my_pca_df2 <- cbind(my_pca_df, coldata)

  # get eigen values (standard deviation squared)
  eigs <- mypca$sdev^2

  # get the percent variance explained by the two PC's
  pc_a_var <- round(((eigs[as.numeric(gsub("PC", "", pc_a))] / sum(eigs)) * 100), 1)
  pc_b_var <- round(((eigs[as.numeric(gsub("PC", "", pc_b))] / sum(eigs)) * 100), 1)

  # construct filename
  filename <- paste0(pc_a, "_", pc_b, "_", "color_by_", color_var, ".", plot_type)

  # plot
  myplot <- ggplot(my_pca_df2, aes_string(pc_a, pc_b, color=color_var)) +
    geom_point(size=5) +
    xlab(paste0(pc_a, ": ", pc_a_var, "% variance")) +
    ylab(paste0(pc_b, ": ", pc_b_var, "% variance")) +
    ggtitle("PCA: Gene Expression") +
    #geom_text(aes_string(label=label_var),hjust=0.7, vjust=-1.1, show.legend = FALSE) +
    geom_text_repel(aes_string(label=label_var)) +
    theme_classic() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  #ggsave(filename=filename)

  return(myplot)
}

#---------------------------------------------------------------------------------------------------------------------------#
