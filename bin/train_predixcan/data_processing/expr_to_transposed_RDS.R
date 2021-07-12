# Reads an expression file in as a dataframe, transposes it, and
# saves it as an RDS object.  If a covariate file is present (as 3rd
# input argument), then it will be used to produce new expression data
# by making a linear model expression ~ covariate, and then pulling the
# residuals as the new expression data.
#
# The input expression file is expected to be tab-delimited, with people
# as columns and genes as row.

# modified by Laura Colbran (1/7/2018) to normalize and split up GTEx expression data
# by tissue, and also to conduct some QC based on read count

library('readr')

# borrowed code from here: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
# tweaked tie handling to "average" to match bioconductor
quantileNormalisation <- function(df){
    df_rank <- apply(df,2,rank,ties.method="average")
    df_sorted <- data.frame(apply(df, 2, sort))
    df_mean <- apply(df_sorted, 1, mean)

    index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
    }

    df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
    rownames(df_final) <- rownames(df)
    return(df_final)
}

# match ids to tissue names based on sample attributes file
idsList <- function(file_path) {
    annotations <- read_delim(gzfile(file_path), "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE,
          comment = "#",col_types = do.call(cols_only, list(SAMPID = col_character(), SMTSD= col_character(), SMAFRZE = col_character())))
    
    tiss_dict <- list()
    tiss_names <- unique(annotations$SMTSD)
    for (i in 1:length(tiss_names)) {
        if (is.na(tiss_names[i])) {
            next
        }
        tiss <- gsub(" - ", "_", tiss_names[i],fixed = TRUE)# remove dashes (these two are to match covariate file name)
        tiss <- gsub(" ", "_", tiss,fixed = TRUE)# remove whitespace
        tiss <- gsub("(", "", tiss,fixed = TRUE)
        tiss <- gsub(")", "", tiss,fixed = TRUE)
        
        ids <- annotations[!is.na(annotations$SMTSD) & annotations$SMTSD == tiss_names[i],]$SAMPID
        ids <- ids[!is.na(ids)]
        for (j in 1:length(ids)) {
            ids[j] <- gsub("-",".",ids[j]) #so that they'll match column names later
        }
        tiss_dict[[tiss]] <- ids
    }
    return(tiss_dict)
}

#filter individual ids and genes by reads
readFilter <- function(tiss_name, ids, read_path) {
    reads <- read.table(read_path, stringsAsFactors = FALSE,
        header = TRUE, row.names = 1,skip=2)
    reads <- reads[intersect(names(reads),ids)]
    
    # summarize across aliquots from same person
    all_ids <- sort(names(reads))
    ids_to_keep <- c(all_ids[1])
    for (i in 2:length(all_ids)) {
        prev <- strsplit(tail(ids_to_keep, n=1),'[.]')[[1]][2]
        curr <- strsplit(tail(all_ids[i], n=1),'[.]')[[1]][2]
        if (prev == curr) {
            if (sum(reads[all_ids[i]]) > sum(reads[tail(ids_to_keep, n=1)])) {
                ids_to_keep[length(ids_to_keep)] <- all_ids[i]
            }
        } else {
            ids_to_keep <- c(ids_to_keep,all_ids[i])
        }
    }
    reads <- reads[ids_to_keep]
    
    # filter for genes with ≥6 reads in at least 10 individuals
    genes <- c()
    for (i in 1:nrow(reads)) {
        read_list <- reads[i,]
        if (length(read_list[which(read_list >= 6)]) >= 10) {
            genes <- c(genes, row.names(reads)[i])
        }
    }    
    return(list(names(reads),genes))
}

makeOutput <- function(tiss_name,targets, rpkm_path, cov_path, out_dir) {
    #print("writing output...")
    rpkm <- read.table(rpkm_path, stringsAsFactors = FALSE,
        header = TRUE, skip=2)
    
    # set rownames to Gene ids, remove non-expression data columns
    rownames(rpkm) <- rpkm[,1]
    rpkm[,-c(1,2)]
    
    # filter to target ids and genes identified based on readFilter    
    rpkm <- rpkm[intersect(names(rpkm),unlist(targets[[1]]))]
    rpkm <- subset(rpkm, row.names(rpkm) %in% targets[[2]])
    
    # filter for genes with >0.1 RPKM in at least 10 individuals
    genes <- c()
    for (i in 1:nrow(rpkm)) {
        rpkm_list <- rpkm[i,]
        if (length(rpkm_list[which(rpkm_list > 0.1)]) >= 10) {
            genes <- c(genes, row.names(rpkm)[i])
        }
    }    
    rpkm <- subset(rpkm, row.names(rpkm) %in% genes)
    #print("Normalizing....")
    # "Expression values were quantile normalized to the average empirical distribution observed across samples."
    # function defined above
    rpkm <- quantileNormalisation(rpkm)
    
    # "For each gene, expression values were inverse quantile normalized to a standard normal distribution across samples."
    # command from eric's link (pg18): https://media.nature.com/original/nature-assets/nature/journal/v490/n7419/extref/nature11401-s1.pdf
    rpkm <- t(rpkm)
    for (i in 1:ncol(rpkm)) {
        norm_dist <- qnorm((rank(rpkm[,i],na.last="keep")-0.5)/sum(!is.na(rpkm[,i])))
        rpkm[,i] <- norm_dist
    }
    
     # truncate ind IDs for matching to genotype
    ids <- row.names(rpkm)
    for (i in 1:length(ids)) {
        ids[i] <- paste(strsplit(ids[i],'[.]')[[1]][1:2],collapse = "-")
    }
    row.names(rpkm) <- ids             
    
    # Correct for covariates if they were provided
    if (!is.na(cov_path)) {
       cov_files <- list.files(cov_path)
       for (i in 1:length(cov_files)) {
           tiss <- strsplit(cov_files[i],"_Analysis")[[1]][1]
           if (tiss_name == tiss) {
              print(tiss_name)
              covariate <- read.table(paste(cov_path,cov_files[i],sep=""), stringsAsFactors = FALSE,
                header = TRUE, row.names = 1)
              covariate <- t(covariate)

              # filter to set of people in Eric's covariate file (accounts for some filtering I couldn't replicate independently)
              ids <- row.names(covariate)
              for (i in 1:length(ids)) {
                  ids[i] <- paste(strsplit(ids[i],'[.]')[[1]][1:2],collapse = "-")
              }
              row.names(covariate) <- ids
              rpkm <- rpkm[rownames(rpkm) %in% rownames(covariate),]

              # confirm same order
              rpkm[ order(row.names(rpkm)),]
              covariate[ order(row.names(covariate)),]

              for (i in 1:length(colnames(rpkm))) {
                fit <- lm(rpkm[,i] ~ covariate) #gene ~ covariates
                rpkm[,i] <- fit$residuals
              }
           }
        }
    }
    saveRDS(rpkm, paste0(paste0(out_dir,tiss_name,""),".RDS",""))
    write.table(t(rpkm), file=paste0(paste0(out_dir,tiss_name,""),".txt",""),quote=F,col.names=T,row.names=T,sep='\t') #for the metadata script
}

argv <- commandArgs(trailingOnly = TRUE)
annotations <- argv[1]
rpkmfile <- argv[2]
readfile <- argv[3]
out_dir <- argv[4]

# handle possible presence of covariate file
covariatepath <- ifelse(length(argv) == 5, argv[5], NA)

tissue_ids <- idsList(annotations)
for (tiss in names(tissue_ids)) {
    print(tiss)
    targets <- readFilter(tiss, tissue_ids[[tiss]],readfile)
    makeOutput(tiss,targets,rpkmfile,covariatepath,out_dir)
}

