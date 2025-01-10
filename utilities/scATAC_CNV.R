#!/usr/bin/env Rscript

library(grid)
library(dplyr)
library(data.table)
library(SummarizedExperiment)
set.seed(1)

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly=TRUE)
FRAG_FILE <- args[[1]]
WINDOW_SIZE <- args[[2]]
STEP_SIZE <- args[[3]]
NEIGHBORS <- args[[4]]
WHITELIST <- args[[5]]
OUTPUT_DIR <- args[[6]]
REFERENCE <- args[[7]]
ROOT_DIR <- args[[8]]

# load in utilities
source(paste0(
    ROOT_DIR, "/utilities/scATAC_CNV_utilities.R"
)

frag_files <- readFragFiles(FRAGS_DIR)
names(frag_files) <- names(frag_files) %>% gsub('_atac.*','',.) %>% gsub('.*/','',.)

whitelist <- read.table(WHITELIST, sep='\t', header=F, comment.char='')[,1]

if (REFERENCE == 'hg19') {
    library(BSgenome.Hsapiens.UCSC.hg19)
    blacklist <- rtracklayer::import.bed(
        paste0(ROOT_DIR, '/references/hg19-blacklist.v2.bed.gz')
    )
    genome <- BSgenome.Hsapiens.UCSC.hg19
} else {
    library(BSgenome.Hsapiens.UCSC.hg38)
    blacklist <- rtracklayer::import.bed(
        paste0(ROOT_DIR, '/references/hg38.blacklist.bed.gz')
    )
    genome <- BSgenome.Hsapiens.UCSC.hg38
}

windows <- makeWindows(genome = genome,
                           blacklist = blacklist %>% subset(width(.) > 1000),
                           windowSize = as.numeric(WINDOW_SIZE),
                           slidingSize = as.numeric(STEP_SIZE))

fragFiles2cnaObj(frag_files = frag_files, 
                     windows = windows,
                     whitelist = whitelist,
                     dir = OUTPUT_DIR,
                     neighbors = as.numeric(NEIGHBORS))

