### qualified variantsubgroup，ACAT-subgroup
library(getopt)
library(magrittr)
library(stringr)
library(readr)
# 
spec <- matrix(
  c("label_file",  "l", 1, "character",  "the file of the label of samples",
    "path_name", "p", 1, "character", "the path prefix for the input/output file",
    "external_MAF", "e", 2, "character", "the external MAF cutoff for lof/mis/syn variant",
    "internal_MAF", "i", 2, "character", "the internal MAF cutoff for lof/mis/syn variant",
    "transcript_type", "t", 1, "character", "which transcript to keep: ccds/refseq",
    "subgroup_file", "s", 1, "character", "the file containing subgroup information",
    "subgroup_col", "n", 1, "character", "the number of column which was used to subgroup",
    "mask_spliceAI", "m", 2, "character", "the mask assigned to a splice variant. mask2:0.8. mask3:0.5. mask4:0.3. mask5:0.1",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)


print(opt)
if( !is.null(opt$help) || is.null(opt$path_name) || is.null(opt$transcript_type) ||
    is.null(opt$label_file) || is.null(opt$subgroup_file) || is.null(opt$subgroup_col) ){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
date()

maxMAF.external <- as.numeric(opt$external_MAF)
maxMAF.internal <- as.numeric(opt$internal_MAF)
transcript.type <- opt$transcript
if(is.null(opt$mask_spliceAI)){
  mask_spliceAI <- "mask2"
}else{
  mask_spliceAI <- opt$mask_spliceAI
}

# subgroup
load(opt$subgroup_file) # subgroup
sub_col <- as.numeric(opt$subgroup_col)
subgroup_name <- colnames(subgroup)[sub_col]
subgroup_sam <- rownames(subgroup)[subgroup[,sub_col]%in%1]
## 
sam.info <- read.table(opt$label_file, sep = "\t", stringsAsFactors = F)
colnames(sam.info) <- c("sam_name", "class")
sam.control <- sam.info$sam_name[sam.info$class=="control"]
sam.all <- c(sam.control, subgroup_sam)
sam.all <- sam.all[sam.all%in%sam.info$sam_name] # 
system(command = paste0("mkdir -p ", dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name)) # 

input.file <- paste0(opt$path_name, ".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".txt")
output.file <- paste0(dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name, "/", basename(opt$path_name),
                      ".",subgroup_name,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".txt")
con <- file(input.file, "r")
ln_num <- 1
ln=readLines(con,n=1)
write_lines(ln, path = output.file)
ln_num <- 2
ln=readLines(con,n=1)
while(length(ln) != 0) {
  ln.ele <- str_split(ln, pattern = "\t") %>% unlist
  ln.sam <- str_split(ln.ele[4], pattern = ";") %>% unlist
  ln.sam.info <- str_split(ln.ele[6], pattern = ";") %>% unlist
  out.sam.idx <- ln.sam %in% sam.all
  out.sam <- ln.sam[out.sam.idx]
  if(length(out.sam)>0){
    out.sam.info <- ln.sam.info[out.sam.idx]
    ln.ele[4] <- paste(out.sam, collapse = ";")
    ln.ele[6] <- paste(out.sam.info, collapse = ";")
    ln.mis.sam <- str_split(ln.ele[7], pattern = ";") %>% unlist
    ln.mis.sam.info <- str_split(ln.ele[8], pattern = ";") %>% unlist
    out.mis.sam.idx <- ln.mis.sam %in% sam.all
    out.mis.sam <- ln.mis.sam[out.mis.sam.idx]
    out.mis.sam.info <- ln.mis.sam.info[out.mis.sam.idx]
    if(length(out.mis.sam)==0){
      out.mis.sam <- ""
      pit.mis.sam.info <- ""
    }
    ln.ele[7] <- paste(out.mis.sam, collapse = ";")
    ln.ele[8] <- paste(out.mis.sam.info, collapse = ";")
    write_lines(paste(ln.ele, collapse = "\t"), path = output.file, append = T)
  }
  ln=readLines(con,n=1)
  ln_num <- ln_num+1
}
close(con)

# tsv
input.file <- paste0(opt$path_name, ".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".tsv")
output.file <- paste0(dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name, "/", basename(opt$path_name),
                      ".",subgroup_name,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".tsv")
print(input.file)
print(output.file)

con <- file(input.file, "r")
ln_num <- 1
ln=readLines(con,n=1)
write_lines(ln, path = output.file)
ln_num <- ln_num+1
ln=readLines(con,n=1)
while(length(ln) != 0) {
  ln.ele <- str_split(ln, pattern = "\t") %>% unlist
  ln.sam <- str_split(ln.ele[4], pattern = ";") %>% unlist
  ln.sam.info <- str_split(ln.ele[6], pattern = ";") %>% unlist
  out.sam.idx <- ln.sam %in% sam.all
  out.sam <- ln.sam[out.sam.idx]
  if(length(out.sam)>0){
    out.sam.info <- ln.sam.info[out.sam.idx]
    ln.ele[4] <- paste(out.sam, collapse = ";")
    ln.ele[6] <- paste(out.sam.info, collapse = ";")
    ln.mis.sam <- str_split(ln.ele[7], pattern = ";") %>% unlist
    ln.mis.sam.info <- str_split(ln.ele[8], pattern = ";") %>% unlist
    out.mis.sam.idx <- ln.mis.sam %in% sam.all
    out.mis.sam <- ln.mis.sam[out.mis.sam.idx]
    out.mis.sam.info <- ln.mis.sam.info[out.mis.sam.idx]
    if(length(out.mis.sam)==0){
      out.mis.sam <- ""
      pit.mis.sam.info <- ""
    }
    ln.ele[7] <- paste(out.mis.sam, collapse = ";")
    ln.ele[8] <- paste(out.mis.sam.info, collapse = ";")
    write_lines(paste(ln.ele, collapse = "\t"), path = output.file, append = T)
  }
  ln=readLines(con,n=1)
  ln_num <- ln_num+1
}
close(con)
date()