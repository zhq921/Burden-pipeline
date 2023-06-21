library(magrittr)
library(stringr)
library(readr)
library(getopt)
# 
spec <- matrix(

  c("input_file",  "i", 1, "character", "relateness2 result file",
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

con <- file(opt$input_file, "r")
ln=readLines(con,n=1)
ln_num <- 1
res <- c()
while(length(ln) != 0){
  if(grepl("-nan$", x = ln)){
    ln_num <- ln_num+1
    ln=readLines(con,n=1)
  }
  ln.info <- str_split(string = ln, pattern = "\t") %>% unlist
  pair1 <- paste(ln.info[1], ln.info[2], sep = ";")
  pair2 <- paste(ln.info[2], ln.info[1], sep = ";")
  if((!pair1 %in% res) & (!pair2 %in% res)){
   res <- c(res, pair1)
  }
  ln_num <- ln_num+1
  ln=readLines(con,n=1)
}

res <- str_split_fixed(res, pattern = ";", n = 2)
res <- unique(res[,1])
res <- cbind(res, "high relatedness")

write.table(res, file = opt$output_file, quote = F, col.names = F, row.names = F, append = T, sep = "\t")