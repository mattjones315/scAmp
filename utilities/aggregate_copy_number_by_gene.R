library(ArchR, quietly=T)
library(ggplot2, quietly=T)
library(patchwork)

# ------------ UTILITIES ------------- #
average_windows_containing_region <- function(cna_data, query_gr, copy_number_array) {
    gr <- rowRanges(cna_data)
    windows <- subsetByOverlaps(gr, query_gr, maxgap=1e5, minoverlap=0L, type='any')

    if (length(windows) == 0) {
        return(NULL)
    }

    cns_subset <- copy_number_array[,windows$name,drop=F]
    return(apply(cns_subset, 1, mean))
}


# ------------- ARGUMENTS -------------- #
args <- commandArgs(trailingOnly=TRUE)
CNA_DIR <- args[[1]]
OUTPUT_DIR <- args[[2]]
REFERENCE <- args[[3]]

# ------------- PIPELINE -------------- #

addArchRGenome(REFERENCE)

if (REFERENCE == 'hg19') {
    genes <- geneAnnoHg19$genes$symbol
} else {
    genes <- geneAnnoHg38$genes$symbol
}

genes <- genes[!is.na(genes)]
arch_object <- loadArchRProject(ARCHR_OBJECT)

files <- list.files(CNA_DIR)
for (file in files) {
    message(paste0("Processing ", file, '...'))

    message("Loading in data...")
    cna_windows <- readRDS(paste0(CNA_DIR, "/", file))

    # creating copy-number array
    cns = assays(cna_windows)$CNs
    rownames(cns) <- rowData(cna_windows)$name
    cns = as.data.frame(t(cns))

    message("Computing copy numbers for gene bodies...")
    cn_list = list()
    x <- 0
    for (gene in genes) {
        if(x %% 1000 == 0){
                message(sprintf("%s of %s", x, length(genes)))
        }
        x <- x + 1

        copy_number <- average_windows_containing_region(cna_windows, subset(genes, symbol == gene), cns)
        if (is.null(copy_number)) {
            next
        }
        
        cn_list[[gene]] <- copy_number
    }

    cn_dataframe <- data.frame(cn_list)
    # names(cn_dataframe) <- genes

    file.stub <- strsplit(file, '.rds', fixed=T)[[1]]
    agg.name <- paste0(OUTPUT_DIR, "/", file.stub, ".all.tsv")

    write.table(cn_dataframe, agg.name, sep='\t')
    message(paste0("Output written to ", agg.name))

}

message("Done.")




