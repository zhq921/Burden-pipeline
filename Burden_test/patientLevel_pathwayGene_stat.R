
library(magrittr)
library(getopt)
library(stringr)
library(data.table)

# 
spec <- matrix(
  c(
    "project_path", "j", 2, "character", "the path for the project",
    "pathway_name", "p", 1, "character", "the pathway to dig into",
    "subgroup_name", "s", 1, "character", "the subgroup name column",
    "variant_type", "v", 1, "character", "the variant type: lof or lof_hc or lof_mis or syn",
    "label_file",  "l", 2, "character",  "the file of the label of samples", 
    "external_MAF", "e", 2, "character", "the external MAF cutoff for lof/mis/syn variant",
    "internal_MAF", "i", 2, "character", "the internal MAF cutoff for lof/mis/syn variant",
    "transcript_type", "t", 2, "character", "which transcript to keep: ccds/refseq",
    "subgroup_file", "f", 2, "character", "the file containing subgroup information",
    "pathway_file", "g", 2, "character", "the RData file for candidate and msigdb pathway info",
    "spl_mask3", "m", 1, "character", "if use spl mask3 data",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
if( !is.null(opt$help) ||
    is.null(opt$variant_type) || is.null(opt$subgroup_name) ){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
print("Start...")

# 
if(is.null(opt$project_path)){ # 
  opt$project_path <- "/share/pub/zhaohq/project/pumch/burden/scoliosis_2021Sep" # CS
}
if(is.null(opt$label_file)){
  opt$label_file <- "./samInfo/sam_class_afFlt.txt"
}
if(is.null(opt$external_MAF)){
  opt$external_MAF <- 1e-04
}
if(is.null(opt$internal_MAF)){
  opt$internal_MAF <- 2.5e-04
}
if(is.null(opt$transcript_type)){
  opt$transcript_type <- "refseq"
}
if(is.null(opt$subgroup_file)){
  opt$subgroup_file <- "./samInfo/subgroup.RData"
}
if(is.null(opt$pathway_file)){
  opt$pathway_file <- "./samInfo/pathway.RData"
}
setwd(opt$project_path) # 
opt$external_MAF <- as.numeric(opt$external_MAF)
opt$internal_MAF <- as.numeric(opt$internal_MAF)
## subgroup
load(opt$subgroup_file) # subgroup
subgroup_name <- colnames(subgroup)[as.numeric(opt$subgroup_name)]
if(!subgroup_name %in% colnames(subgroup)){
  # ，
  paste("Error: wrong subgroup name!", opt$subgroup_name, "does not exist!") %>% print()
  q()
}
subgroup_sam <- rownames(subgroup)[subgroup[,subgroup_name]==1]

## pathway
load(opt$pathway_file)
pathway_name <- opt$pathway_name
if(!pathway_name %in% names(pathway)){
  # ，
  paste("Error: wrong pathway name!", pathway_name, "does not exist!") %>% print()
  q()
}
pathway_genes <- pathway[[pathway_name]]
if(is.null(pathway_genes)){print("No genes in this pathway")}
pathway_name <- gsub("[ ]+", "_", pathway_name)

# input and output file
if(opt$spl_mask3 == "T"){
  input.file <- paste0("./burdenTest/subgroup/qualified_variant/", subgroup_name, "/", basename(opt$project_path),
                       ".indelRmrptSNPraw.",subgroup_name,".",opt$transcript_type,".e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.txt")
}else{
  input.file <- paste0("./burdenTest/subgroup/qualified_variant/", subgroup_name, "/", basename(opt$project_path),
                       ".indelRmrptSNPraw.",subgroup_name,".",opt$transcript_type,".e",opt$external_MAF,".i",opt$internal_MAF,".txt")
  
}
if(!opt$variant_type %in% c("lof_mis", "lof_hc", "lof", "syn")){

  paste("Error: wrong variant type!", opt$variant_type, "does not exist! Only lof,lof_hc,lof_mis,syn are allowed!") %>% print()
  q()
}
system(command = paste0("mkdir -p ", "./burdenTest/pathway_specified/", 
                        subgroup_name, "/", "e",opt$external_MAF,
                        ".i",opt$internal_MAF, "/", opt$variant_type )) # 
if(opt$variant_type=="lof_mis"){ # lof_mis，lof+misburden
  output.file <- paste0("./burdenTest/pathway_specified/", subgroup_name, "/", 
                        "e",opt$external_MAF, ".i",opt$internal_MAF, "/", opt$variant_type, "/", pathway_name, ".",
                        opt$transcript_type,".lof_mis.e",opt$external_MAF,".i",opt$internal_MAF,".sample.tsv")
  masks <- c("mask1","mask2","mask3","mask4","mask5",
             "mask1m","mask2m","mask3m","mask4m","mask5m")
  mask_weight <- c(1, 0.8, 0.5, 0.3, 0.1, 
                   1, 0.8, 0.5, 0.3, 0.1)
}else{
  if(opt$variant_type=="lof"){
    output.file <- paste0("./burdenTest/pathway_specified/", subgroup_name, "/", 
                          "e",opt$external_MAF, ".i",opt$internal_MAF, "/", opt$variant_type, "/", pathway_name, ".",
                          opt$transcript_type,".lof.e",opt$external_MAF,".i",opt$internal_MAF,".sample.tsv")
    masks <- c("mask1","mask2")
    mask_weight <- c(1, 0.8)
  }else{
    if(opt$variant_type=="syn"){
      output.file <- paste0( "./burdenTest/pathway_specified/", subgroup_name, "/", 
                            "e",opt$external_MAF, ".i",opt$internal_MAF, "/", opt$variant_type, "/", pathway_name, ".",
                            opt$transcript_type,".syn.e",opt$external_MAF,".i",opt$internal_MAF,".sample.tsv")
      masks <- c("syn")
      mask_weight <- c(1)
    }else{
      if(opt$variant_type=="lof_hc"){
        output.file <- paste0("./burdenTest/pathway_specified/", subgroup_name, "/", 
                              "e",opt$external_MAF, ".i",opt$internal_MAF, "/", opt$variant_type, "/", pathway_name, ".",
                              opt$transcript_type,".lof_hc.e",opt$external_MAF,".i",opt$internal_MAF,".sample.tsv")
        masks <- c("mask1")
        mask_weight <- c(1)
      }
    }
  }
}

print(opt)
print("Input file:")
print(input.file)
print("Output file:")
print(output.file)

## 
sam.info <- read.table(opt$label_file, sep = "\t", stringsAsFactors = F)
colnames(sam.info) <- c("sam_name", "class")
sam.control <- sam.info$sam_name[sam.info$class=="control"]
sam.case <- sam.info$sam_name[sam.info$class=="case"]
# print(paste0("the number of subgroup sample: ", length(subgroup_sam)))
# print(paste0("the number of case sample: ", length(sam.case)))
sam.case <- intersect(sam.case, subgroup_sam)
print(paste0("the number of subgroup & case sample: ", length(sam.case)))
print(paste0("the number of control sample: ", length(sam.control)))
sam.n.ctrl <- length(sam.control)
sam.n.case <- length(sam.case)
sam.info <- data.frame(sam_name = c(sam.case, sam.control),
                       class = c(rep("case", sam.n.case), rep("control", sam.n.ctrl)),
                       stringsAsFactors = F
)
case.idx <- sam.info$class=="case"
cont.idx <- sam.info$class=="control"

# ，
vrt2sam <- fread(input.file, sep = "\t", stringsAsFactors = F, header = T, data.table = F) 
vrt2sam <- vrt2sam[vrt2sam[,2]%in%masks,] # masks

## ，，hash
sam.info <- data.frame(sam.info, gene_num = 0, wtd_gene_num = 0, 
                       genes = NA, masks = NA, stringsAsFactors = F)
patient.info <- list()
for(i in 1:nrow(vrt2sam)){ # 
  ln.ele  <- vrt2sam[i,] %>% unlist()
  ln.mask <- ln.ele[2]
  ln.gene <- str_split(ln.ele[3], pattern = ";") %>% unlist # 
  
  # pathway
  ln.gene <- ln.gene[ln.gene %in% pathway_genes]
  if(length(ln.gene) == 0){next} # pathway，
  
  ln.mtt.sam <- str_split(ln.ele[4], pattern = ";") %>% unlist # 
  for(j in 1:length(ln.gene)){ # 
    gene <- ln.gene[j]
    for(k in 1:length(ln.mtt.sam)){ # 
      pt <- ln.mtt.sam[k]
      if(!pt %in% names(patient.info)){ # 
        patient.info <- c(patient.info, list(matrix(c(gene, ln.mask), nrow = 1)))
        names(patient.info)[length(patient.info)] <- pt
      }else{
        if(!gene %in% patient.info[[pt]][,1]){ # 
          patient.info[[pt]] <- rbind(patient.info[[pt]], c(gene, ln.mask))
        }else{ # ，mask
          if(ln.mask < patient.info[[pt]][match(gene, patient.info[[pt]][,1]),2]){ # mask
            patient.info[[pt]][match(gene, patient.info[[pt]][,1]),2] <- ln.mask
          }
        }
      }
    }
  }
}

# 
for(i in 1:length(patient.info)){
  pt <- names(patient.info)[i]
  pt.info <- patient.info[[pt]]
  pt.genes <- paste(pt.info[,1], collapse = ";")
  pt.masks <- paste(pt.info[,2], collapse = ";")
  pt.gene.num <- nrow(pt.info)
  pt.wtd.gene.num <- mask_weight[match(pt.info[,2], masks)] %>% sum()
  sam.info[match(pt, sam.info$sam_name), 3:6] <- c(pt.gene.num, pt.wtd.gene.num, pt.genes, pt.masks)
}
sam.info$gene_num <- sam.info$gene_num %>% as.numeric()
sam.info$wtd_gene_num <- sam.info$wtd_gene_num %>% as.numeric()

#  
sam.info.stat <- table(sam.info$class, sam.info$gene_num)
err <- class(try(fisher.test(sam.info.stat, workspace = 1e07), silent = T)) 
if("try-error" %in% err){ # fisher
  sam.info.stat.p <- (sam.info.stat %>% 
                        fisher.test(., simulate.p.value = T, B = 1000000))$p.value %>% 
    signif(., digits = 3)
}else{
  sam.info.stat.p <- (sam.info.stat %>% fisher.test(., workspace = 1e07))$p.value %>% signif(., digits = 3)
}
sam.info.stat.perc <- sam.info.stat *100 / c(sam.n.case, sam.n.ctrl)
sam.info.stat.perc <- apply(sam.info.stat.perc, 2, function(x) {round(x, digits = 2)})
#  
sam.info.wtd.stat <- table(sam.info$class, sam.info$wtd_gene_num %>% ceiling(.) %>% paste0("<",.))
err <- class(try(fisher.test(sam.info.wtd.stat, workspace = 1e07), silent = T)) 

if("try-error" %in% err){ # fisher
  sam.info.wtd.stat.p <- (sam.info.wtd.stat %>% 
                            fisher.test(., simulate.p.value = T, B = 1000000))$p.value %>% 
    signif(., digits = 3)
}else{
  sam.info.wtd.stat.p <- (sam.info.wtd.stat %>% 
                            fisher.test(., workspace = 1e07))$p.value %>% 
    signif(., digits = 3)
}
sam.info.wtd.stat.perc <- sam.info.wtd.stat *100 / c(sam.n.case, sam.n.ctrl)
sam.info.wtd.stat.perc <- apply(sam.info.wtd.stat.perc, 2, function(x) {round(x, digits = 2)})

# warning，append header
options(warn = -1)

# 
write.table("Distribution by gene counts:", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F)
write.table(sam.info.stat, file = output.file, sep = "\t", quote = F, col.names = NA, row.names = T, append = T)
write.table("Distribution by gene counts (%):", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table(sam.info.stat.perc, file = output.file, sep = "\t", quote = F, col.names = NA, row.names = T, append = T)
write.table(paste0("P-value: ", sam.info.stat.p), file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table("", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table("Distribution by weighted gene counts:", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table(sam.info.wtd.stat, file = output.file, sep = "\t", quote = F, col.names = NA, row.names = T, append = T)
write.table("Distribution by weighted gene counts (%):", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table(sam.info.wtd.stat.perc, file = output.file, sep = "\t", quote = F, col.names = NA, row.names = T, append = T)
write.table(paste0("P-value: ", sam.info.wtd.stat.p), file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table("", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table("Detail info:", file = output.file, sep = "\t", quote = F, col.names = F, row.names = F, append = T)
write.table(sam.info, file = output.file, sep = "\t", quote = F, col.names = T, row.names = F, append = T)
print("End!")



