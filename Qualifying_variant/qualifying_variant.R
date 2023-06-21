library(magrittr)
library(stringr)
library(readr)
library(getopt)
source("/mnt/work/research/zhaohq/burden/code/function.R")
print(getwd())

spec <- matrix(
  c("input_vcf_file",  "v", 1, "character", "vcf file after vep",
    "output_file", "o", 1, "character",  "output_file_prefix",
    "external_lof_maf",  "e", 1, "double",  "max external MAF for lof",
    "internal_lof_maf",  "i", 1, "double",  "max internal MAF for lof",
    "transcript", "t", 1, "character", "which transcript to keep: ccds/refseq",
    "mask_spliceAI", "m", 2, "character", "the mask assigned to a splice variant. mask2:0.8. mask3:0.5. mask4:0.3. mask5:0.1",
    "help",   "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
print(opt)
if( !is.null(opt$help) || is.null(opt$input_vcf_file) || is.null(opt$output_file)  || 
    is.null(opt$external_lof_maf) || is.null(opt$internal_lof_maf) || is.null(opt$transcript) ){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
date()

con <- file(opt$input_vcf_file, "r")
input.nrow <- system(command = paste("wc -l",opt$input_vcf_file), intern = T) %>% 
  strsplit(., split = " ") %>% unlist() %>% extract(.,1) %>% 
  as.numeric() # ，
out.file <- opt$output_file
maxMAF.external <- as.numeric(opt$external_lof_maf)
maxMAF.internal <- as.numeric(opt$internal_lof_maf)
transcript.type <- opt$transcript
if(is.null(opt$mask_spliceAI)){
  mask_spliceAI <- "mask2"
}else{
  mask_spliceAI <- opt$mask_spliceAI
}

# 
ln_num <- 0
while(!exists("vcf.col.nam")) {
  ln_num <- ln_num+1
  ln=readLines(con,n=1)
  if(str_detect(ln, "^##INFO=<ID=CSQ")){
    vep.col.nam <- ln %>% str_split(., pattern = ": ") %>% unlist %>% extract(.,2) %>%
      str_split(., pattern = "\"") %>% unlist %>% extract(.,1) %>%
      str_split(., pattern = "[|]") %>% unlist
    vep.col.len <- length(vep.col.nam)
    database.maxaf.idx <- which(vep.col.nam=="MAX_AF") # highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD
    gnomad.af.idx <- which(vep.col.nam=="gnomAD_exomes_controls_POPMAX_AF") # Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry) in the controls subset
    gnomadg.af.idx <- which(vep.col.nam=="gnomAD_genomes_POPMAX_AF")
    canonical.idx <- which(vep.col.nam=="CANONICAL")
    ccds.idx <- which(vep.col.nam=="CCDS")
    symbol.idx <- which(vep.col.nam=="SYMBOL")
    vrt.csq.idx <- which(vep.col.nam=="Consequence")
    spliceai.idx <- str_which(vep.col.nam, pattern = "SpliceAI_pred_DS")
    spliceai.symbol.idx <- which(vep.col.nam=="SpliceAI_pred_SYMBOL")
    loftee.idx <- which(vep.col.nam=="LoF")
    cadd.idx <- which(vep.col.nam=="CADD_PHRED")
    vest.idx <- which(vep.col.nam=="VEST4_rankscore")
    dann.idx <- which(vep.col.nam=="DANN_rankscore")
    eRaw.idx <- which(vep.col.nam=="Eigen-raw_coding_rankscore")
    ePCraw.idx <- which(vep.col.nam=="Eigen-PC-raw_coding_rankscore")
    fathmm.idx <- which(vep.col.nam=="FATHMM_pred")
    fathmm.mkl.idx <- which(vep.col.nam=="fathmm-MKL_coding_pred")
    provean.idx <- which(vep.col.nam=="PROVEAN_pred")
    metasvm.idx <- which(vep.col.nam=="MetaSVM_pred")
    metalr.idx <- which(vep.col.nam=="MetaLR_pred")
    mcap.idx <- which(vep.col.nam=="M-CAP_score")
    primate.idx <- which(vep.col.nam=="PrimateAI_pred")
    revel.idx <- which(vep.col.nam=="REVEL_score")
    # ，
    write_lines(str_c("variant name", "min mask level", "gene symbol", "mutated sample", "genotype format", "mutated sample genotype", "missing sample", "missing sample genotype", sep = "\t"), 
                file = paste0(out.file,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".txt"), append = F)
    write_lines(str_c("variant name", "mask level", "gene symbol", "mutated sample", "genotype format", "mutated sample genotype", "missing sample", "missing sample genotype", "variant info", str_c(vep.col.nam, collapse = "\t"), sep = "\t"), 
                file = paste0(out.file,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".tsv"), append = F)
  }  # VEP,vep
  
  if(str_detect(ln, "^#CHROM")){
    vcf.col.nam <- ln %>% str_split(., pattern = "#") %>% unlist %>% extract(.,2) %>%
      str_split(., pattern = "\t") %>% unlist
    info.idx <- which(vcf.col.nam=="INFO")
    format.idx <- which(vcf.col.nam=="FORMAT")
  } # vcf
}

# qualifying variant
while(ln_num < input.nrow){
  ln_num <- ln_num+1
  ln=readLines(con,n=1)
  vcf.col <- ln %>% str_split(., pattern = "\t") %>% unlist # vcf
  info.col <- vcf.col[info.idx] %>% str_split(., pattern = ";") %>% unlist # vcfINFO
  af.internal <- info.col %>% extract(., str_which(info.col, pattern = "^AF=")) %>% 
    str_split(., pattern = "=") %>% unlist %>% extract(., 2) %>% 
    str_split(., pattern = ",") %>% unlist %>% extract(., 1) %>% as.numeric()
  if((af.internal >= maxMAF.internal) | (af.internal==0)){
    next
  } # in house；AB，missing，，
  vrt.vep <- info.col %>% extract(., str_which(info.col, pattern = "^CSQ=")) %>% 
    str_replace(., pattern = "CSQ=", replacement = "") %>% 
    str_split(., pattern = ",") %>% unlist %>% 
    str_split_fixed(., pattern = "[|]", n = vep.col.len) # vep
  if(nrow(vrt.vep)==0){next} # multiallelic，G>*，vep
  # 
  flt <- get_max_af(vrt.vep[, database.maxaf.idx]) < maxMAF.external
  if(!flt){next} # external
  vrt.info.withoutvep <- info.col %>% extract(., str_which(info.col, pattern = "^CSQ=", negate = T)) %>%
    str_c(., collapse = "|") # vep，|，
  # transcript selection
  if(transcript.type=="ccds"){
    vrt.vep <- vrt.vep[str_detect(vrt.vep[,ccds.idx], pattern = "CCDS"),,drop=F] # CCDS transprict
  }else{
    if(transcript.type=="refseq"){
      vrt.vep <- vrt.vep[str_which(vrt.vep[,canonical.idx], pattern = "YES"),,drop=F] # canonical transprict
    }
  }
  if(nrow(vrt.vep)==0){next} # ccds
  # 
  vrt.type <- apply(vrt.vep, 1, function(x) get_variant_type_revel(dat = x, spliceMask = mask_spliceAI))
  if(str_detect(vrt.type, pattern = "t-other") %>% all()){next}# other，
  vrt.vep.selected <- apply(vrt.vep[vrt.type!="t-other",,drop = F], 1, function(x) {str_c(x, collapse = "\t")}) # 
  vrt.type.selected <- vrt.type[vrt.type!="t-other"]
  variant <- str_c(vcf.col[c(1,2,4,5)], collapse = "_")
  mtt.sam.nam <- vcf.col.nam[str_which(vcf.col, pattern = "^1/1|^0/1|^1/0")] %>% str_c(., collapse = ";") # missing，gwasSKATimpute
  mtt.sam.info <- vcf.col[str_which(vcf.col, pattern = "^1/1|^0/1|^1/0")] %>% str_c(., collapse = ";") # 1/01/2multiallelic，，pipeline，1/0，
  mis.sam.nam <- vcf.col.nam[str_which(vcf.col, pattern = "^[.]/[.]")] %>% str_c(., collapse = ";") # missing，gwasSKATimpute
  mis.sam.info <- vcf.col[str_which(vcf.col, pattern = "^[.]/[.]")] %>% str_c(., collapse = ";") # 1/01/2multiallelic，，pipeline，1/0，
  if(length(mis.sam.nam)==0){mis.sam.nam <- ""; mis.sam.info <- "";}
  # 
  # vrt.type.idx <- vrt.type %in% c("mask1","mask2","mask3","mask4","mask5",
  #                                 "mask1m","mask2m","mask3m","mask4m","mask5m", "syn")
  # 20220820，mask6n，revel
  vrt.type.idx <- vrt.type %in% c("mask1","mask2","mask3","mask4","mask5",
                                  "mask1m","mask2m","mask3m","mask4m","mask5m", 
                                  "mask6n","syn")
  if(any(vrt.type.idx)){
    vrt.symbol.all <- vrt.vep[,symbol.idx][vrt.type!="t-other"]
    vrt.symbol.unique <- vrt.symbol.all %>% unique()
    vrt.type.out <- min(vrt.type)
    # ，mask，
    if(length(vrt.symbol.unique)>1){ # 
      vrt.type.out <- c()
      for(s in 1:length(vrt.symbol.unique)){
        vrt.type.out <- c(vrt.type.out, 
                          min(vrt.type[vrt.vep[,symbol.idx]==vrt.symbol.unique[s]]))
      }
    }
    # ，fisher,，
    write_lines(str_c(variant, vrt.type.out, vrt.symbol.unique, mtt.sam.nam, vcf.col[9], mtt.sam.info, mis.sam.nam, mis.sam.info, sep = "\t"), 
                path = paste0(out.file,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".txt"), append = T)
    # ，tab，，
    write_lines(str_c(variant, vrt.type.selected, vrt.symbol.all, mtt.sam.nam, vcf.col[9], mtt.sam.info, mis.sam.nam, mis.sam.info, vrt.info.withoutvep, vrt.vep.selected, sep = "\t"), 
                path = paste0(out.file,".",transcript.type,".e",maxMAF.external,".i",maxMAF.internal,".spl",mask_spliceAI,".tsv"), append = T)
    next
  }
}
close(con)
date()





