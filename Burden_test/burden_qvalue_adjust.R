library(magrittr) # for %>%
library(getopt)
library(data.table)
library(qvalue)

# 
spec <- matrix(
  c("path_name", "p", 1, "character", "the path for the input files.",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$path_name) ){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}

files <- list.files(opt$path_name, pattern = "*.perm$", recursive = T, full.names = T)
for(i in 1:length(files)){
 dat <- fread(input = files[i], header = T, data.table = F)
 tt <- try(burden.p1.qvalue <- qvalue(dat$burden.p1)$qvalues)
 if("try-error" %in% class(tt)){
   burden.p1.qvalue <- qvalue(dat$burden.p1, pi0 = 1)$qvalues
 }
 tt <- try(skat.p1.qvalue <- qvalue(dat$skat.p1)$qvalues)
 if("try-error" %in% class(tt)){
   skat.p1.qvalue <- qvalue(dat$skat.p1, pi0 = 1)$qvalues
 }
 tt <- try(acat.p1.qvalue <- qvalue(dat$acatv.p1)$qvalues)
 if("try-error" %in% class(tt)){
   acat.p1.qvalue <- qvalue(dat$acatv.p1, pi0 = 1)$qvalues
 }
 dat <- cbind(dat, burden.p1.qvalue = burden.p1.qvalue, skat.p1.qvalue = skat.p1.qvalue,
              acat.p1.qvalue = acat.p1.qvalue)
 write.table(dat, 
             file = paste0(files[i], ".qvalue"),
             quote = F, col.names = T, row.names = F,
             sep = "\t")
}