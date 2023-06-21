
library(stringr)
library(magrittr)
library(readr)
library(getopt)


spec <- matrix(

  c("input_vcf_file",  "v", 1, "character", "vcf file after vep",
    "output_file", "o", 1, "character",  "output_file_prefix",
    "minAB",  "i", 1, "double",  "min AB",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
print(opt)
if( !is.null(opt$help) || is.null(opt$input_vcf_file) || is.null(opt$output_file)  || 
    is.null(opt$minAB)){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

date()
con <- file(opt$input_vcf_file, "r")
out.file <- opt$output_file
# maxAB <- 1
minAB <- opt$minAB

ln_num <- 1
ln=readLines(con,n=1)
while(length(ln) != 0) {
  if(exists("vcf.col.nam")){
    vcf.col <- ln %>% str_split(., pattern = "\t") %>% unlist # vcf
    
    format.col <- vcf.col[format.idx] %>% str_split(., pattern = ":") %>% unlist # vcfFORMAT
    ad.idx <- str_which(format.col, pattern = "AD")
    mtt.idx <- grep("^0/1|^1/0", vcf.col)
    vcf.col[mtt.idx] <- 
      sapply(mtt.idx, simplify = T, FUN = function(x){
        dat <- vcf.col[x]
        ad <- str_split(dat, pattern = ":", simplify = T) %>% 
          extract(., 1, ad.idx) %>% 
          str_split(., pattern = ",") %>% unlist %>%
          as.numeric()
        ab <- ad[2]/sum(ad)
        if((ab < minAB) | is.na(ab)){ # AB0.25，./. ；0/1:.:0,0，
          str_sub(dat, 1, 3) <- "./."
        }
        return(dat)
      })
    ln.new <- str_c(vcf.col, collapse = "\t")
    write_lines(ln.new, path = out.file, append = T)
  }else{
    write_lines(ln, path = out.file, append = T)
    if(str_detect(ln, "^#CHROM")){
      vcf.col.nam <- ln %>% str_split(., pattern = "#") %>% unlist %>% extract(.,2) %>%
        str_split(., pattern = "\t") %>% unlist
      info.idx <- which(vcf.col.nam=="INFO")
      format.idx <- which(vcf.col.nam=="FORMAT")
      sam.stt.idx <- format.idx + 1
      sam.end.idx <- length(vcf.col.nam)
    } # vcf
  }
  ln=readLines(con,n=1)
  ln_num <- ln_num+1
}
close(con)
date()