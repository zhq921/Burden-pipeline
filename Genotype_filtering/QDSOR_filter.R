# QD/SOR filter qd2，sor9

library(stringr)
library(magrittr)
library(readr)
library(getopt)

spec <- matrix(

  c("input_vcf_file",  "v", 1, "character", "vcf file after vep",
    "output_file", "o", 1, "character",  "output_file_prefix",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
print(opt)
if( !is.null(opt$help) || is.null(opt$input_vcf_file) || is.null(opt$output_file)){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

date()
con <- file(opt$input_vcf_file, "r")
out.file <- opt$output_file

ln_num <- 1
ln=readLines(con,n=1)
while(length(ln) != 0) {
  if(exists("vcf.col.nam")){
    vcf.col <- ln %>% str_split(., pattern = "\t") %>% unlist # vcf
    info.col <- vcf.col[info.idx] %>% str_split(., pattern = ";") %>% unlist # vcfINFO
    # indel 
    if((nchar(vcf.col[ref.idx]) > 10) |
       (nchar(vcf.col[alt.idx]) > 10)
       ){
      ln=readLines(con,n=1)
      ln_num <- ln_num+1
      next()} # ref_length < 10 & alt_length < 10
    
    # QD, SOR
    qd <- info.col %>% extract(., str_which(info.col, pattern = "^QD=")) %>% 
      str_split(., pattern = "=") %>% unlist %>% extract(., 2) %>% 
      str_split(., pattern = ",") %>% unlist %>% extract(., 1) %>% as.numeric()
    sor <- info.col %>% extract(., str_which(info.col, pattern = "^SOR=")) %>% 
      str_split(., pattern = "=") %>% unlist %>% extract(., 2) %>% 
      str_split(., pattern = ",") %>% unlist %>% extract(., 1) %>% as.numeric()
    if(length(qd)==0 | length(sor)==0){
      print(paste0("line number: ", ln_num, " - without QD/SOR"))
      print(vcf.col[1:9])
      ln=readLines(con,n=1)
      ln_num <- ln_num+1
      next()
    }
    if(qd <= 2 | sor >= 9){ # qd2，sor9
      ln=readLines(con,n=1)
      ln_num <- ln_num+1
      next()
    }
    write_lines(ln, path = out.file, append = T)
  }else{
    write_lines(ln, path = out.file, append = T)
    if(str_detect(ln, "^#CHROM")){
      vcf.col.nam <- ln %>% str_split(., pattern = "#") %>% unlist %>% extract(.,2) %>%
        str_split(., pattern = "\t") %>% unlist
      ref.idx <- which(vcf.col.nam=="REF")
      alt.idx <- which(vcf.col.nam=="ALT")
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