# reads normalized expression bed files downloaded from GTEx portal, reformats,
# then saves as RDS

library('readr')

makeOutput <- function(tiss_name, norm_path, out_dir) {
    #print("writing output...")
    norm <- read.table(norm_path, stringsAsFactors = FALSE,
        header = TRUE,comment.char="z",check.names=FALSE)

    # set rownames to Gene ids, remove non-expression data columns
    rownames(norm) <- norm[,1]
    norm <- norm[,-c(1)]
    head(norm)

    saveRDS(t(norm), paste0(paste0(out_dir,tiss_name,""),".RDS",""))
    write.table(norm, file=paste0(paste0(out_dir,tiss_name,""),".txt",""),quote=F,col.names=T,row.names=T,sep='\t') #for the metadata script
}

argv <- commandArgs(trailingOnly = TRUE)
tiss<- argv[1]
norm_file <- paste(argv[2],"v8.exp.residual.gz",sep=".")
out_dir <- argv[3]

print(tiss)
makeOutput(tiss,norm_file,out_dir)
