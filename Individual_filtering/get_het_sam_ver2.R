
library(getopt)

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

dat <- read.table(file = opt$input_file, header = T, row.names = 1)
dat.mean <- mean(dat$F)
dat.sd <- sd(dat$F)
res.lowhet <- cbind(rownames(dat)[dat$F>(dat.mean+3*dat.sd)], "low heterozygosity")
res.highhet <- cbind(rownames(dat)[dat$F<(dat.mean-3*dat.sd)], "high heterozygosity")
res <- res.lowhet # rbind(res.lowhet, res.highhet)
write.table(res, file = opt$output_file, 
            row.names = F, quote = F, col.names = F, 
            append = T, sep = "\t")
