#!/usr/bin/Rscript

# Take multiple ReadsPerGene.tab files from STAR and combine them into one table
# We will select column 4 (stranded reverse) or column 2 (unstranded) read counts

library(parallel)
library(optparse)
library(readxl)

# Command Line Arguments
# -d : path to directory of STAR's .ReadsPerGene files
# -o : desired outfile (counts_table.txt)
# -s : strand info, which column of STAR ReadsPerGene files to keep (usually 2 or 4)
# -m : path to metadata excel file
#      first column must be desired sample name
#       must contain "fastq" column with fastq file prefix to remove/match
       
option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="path to ReadsPerGene files", metavar="character"),
  make_option(c("-s", "--strand"), type="numeric", default=NULL,
              help="number signifying which column of counts to keep", metavar="numeric"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="path to metadata excel file, first column is desired sample name", metavar="character"),	      
  make_option(c("-o", "--outfile"), type="character", default="counts_table.txt",
              help="desired outfile (counts_table.txt)", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Input Args Checks
if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("Missing input file directory path.n", call.=FALSE)
}
if (is.null(opt$strand)){
  print_help(opt_parser)
  stop("Missing strand column.n", call.=FALSE)
}
#if (is.null(opt$metadata)){
#  print_help(opt_parser)
#  stop("Missing metadata file.n", call.=FALSE)
#}

# create a list that contains the path to all input files
counts_list <- list.files(opt$directory, 'ReadsPerGene.out.tab$', full=T)
names(counts_list) <- counts_list
names(counts_list) <- gsub('.star.ReadsPerGene.out.tab$','',names(counts_list))
names(counts_list) <- gsub(opt$directory, '',names(counts_list))
names(counts_list) <- gsub("/", '',names(counts_list))

# replace the current names (the fastq prefix) with desired sample names from metadata excel file
# read in the metadata excel file
#my_meta <- as.data.frame(read_excel(opt$metadata))
#colnames(my_meta)[1] <- "sample"
# create vector of sample names from matching fastq with sample in metadata file
#my_samples <- as.character(sapply(names(counts_list), function(x) my_meta[my_meta$fastq == x, "sample"])

#names(counts_list) <- my_meta$sample

# read in each counts file as a data frame and store in a list of df's
counts_dfs <- mclapply(counts_list, function(x) read.table(x), mc.cores=8)


# convert the correct counts column to numeric
for (i in 1:length(counts_dfs)) {counts_dfs[[i]][,opt$strand] <- as.numeric(counts_dfs[[i]][,opt$strand])}

# change the name of the correct column to be the sample name
for (i in 1:length(counts_dfs)) {colnames(counts_dfs[[i]])[opt$strand] <- names(counts_dfs)[i]}

# change the name of column 1 to "GeneID"
for (i in 1:length(counts_dfs)) {colnames(counts_dfs[[i]])[1] <- "GeneID"}

# drop unnecessary counts columns
for (i in 1:length(counts_dfs)) {counts_dfs[[i]] <- subset(counts_dfs[[i]], select = colnames(counts_dfs[[i]][c(1,opt$strand)]))}

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
my_table2 <- my_table[ ,match(names(counts_list), colnames(my_table))]
my_table2$GeneID <- rownames(my_table2)
my_table2 <- my_table2[ ,c(ncol(my_table2),1:ncol(my_table2)-1)]

# export the final table
write.table(my_table2, opt$outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
