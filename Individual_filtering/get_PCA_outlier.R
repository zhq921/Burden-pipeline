
library(getopt)
# 
spec <- matrix(

  c("input_file",  "i", 1, "character", "heterozygosity result file",
    "output_file", "o", 1, "character",  "output file",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
print(opt)
if(!is.null(opt$help) || 
   is.null(opt$input_file) || 
   is.null(opt$output_file)
){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

get_outlier <- function(value){
  value.mean <- mean(value)
  value.sd <- sd(value)
  res <- which(value < value.mean-3*value.sd)
  res <- c(res, which(value > value.mean+3*value.sd))
  res <- sort(unique(res))
  return(res)
}

dat <- read.table(file = opt$input_file, row.names = 1)
dat <- dat[,-1]
res <- rownames(dat)[unique(c(get_outlier(dat[,1]), get_outlier(dat[,2])))]
res <- cbind(res, "PCA outlier")
write.table(res, file = opt$output_file, 
            row.names = F, quote = F, col.names = F, 
            append = T, sep = "\t")

