

library(magrittr)
library(hash)
library(getopt)
library(stringr)

# 
spec <- matrix(
  c("label_file",  "l", 1, "character",  "the file of the label of samples", 
    "path_name", "p", 1, "character", "the path for the input file",
    "external_MAF", "e", 2, "character", "the external MAF cutoff for lof/mis/syn variant",
    "internal_MAF", "i", 2, "character", "the internal MAF cutoff for lof/mis/syn variant",
    "variant_type", "v", 2, "character", "the variant type:lof or lof_hc or lof_mis or syn",
    "covariate_PCA", "a", 2, "character", "whether to include covariate using PCA result (PC1:3) with 1000G: the PCA result file path",
    "transcript_type", "t", 1, "character", "which transcript to keep: ccds/refseq",
    "subgroup_file", "s", 1, "character", "the file containing subgroup information",
    "subgroup_col", "n", 1, "character", "the number of column which was used to subgroup",
    "pathway_file", "g", 1, "character", "the RData file for candidate and msigdb pathway info",
    "splmask3", "m", 2, "character", "if use splmask3 data",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$path_name) || is.null(opt$transcript_type) ||
    is.null(opt$label_file) || is.null(opt$pathway_file) ||
    is.null(opt$subgroup_file) || is.null(opt$subgroup_col) ){
  # ... 
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
date()

# subgroup
load(opt$subgroup_file) # subgroup
sub_col <- as.numeric(opt$subgroup_col)
subgroup_name <- colnames(subgroup)[sub_col]
subgroup_sam <- rownames(subgroup)[subgroup[,sub_col]==1]
system(command = paste0("mkdir -p ", dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/",subgroup_name)) # 

opt$external_MAF <- as.numeric(opt$external_MAF)
opt$internal_MAF <- as.numeric(opt$internal_MAF)

if(is.null(opt$splmask3)){
  opt$splmask3 <- "mask2"
}else{
  opt$splmask3 <- opt$splmask3
}

if(opt$splmask3 == "splmask3"){
  input.file <- paste0(dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name, "/", basename(opt$path_name),
                       ".",subgroup_name,".",opt$transcript_type,".e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.txt")
  
}else{
  input.file <- paste0(dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name, "/", basename(opt$path_name),
                       ".",subgroup_name,".",opt$transcript_type,".e",opt$external_MAF,".i",opt$internal_MAF,".splmask2.txt")
  
}
if(opt$variant_type=="lof_mis"){ # lof_mis，lof+misburden
  if(opt$splmask3 == "splmask3"){
    output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                          ".",subgroup_name,".", opt$transcript_type,".lof_mis.e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.pathway.logistic")
  }else{
    output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                          ".",subgroup_name,".", opt$transcript_type,".lof_mis.e",opt$external_MAF,".i",opt$internal_MAF,".splmask2.pathway.logistic")
    
  }
  masks <- c("mask1","mask2","mask3","mask4","mask5",
             "mask1m","mask2m","mask3m","mask4m","mask5m")
  mask_weight <- c(1, 0.8, 0.5, 0.3, 0.1, 
                   1, 0.8, 0.5, 0.3, 0.1)
}else{
  if(opt$variant_type=="lof"){
    if(opt$splmask3 == "splmask3"){
      output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                            ".",subgroup_name, ".", opt$transcript_type,".lof.e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.pathway.logistic")
      
    }else{
      output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                            ".",subgroup_name, ".", opt$transcript_type,".lof.e",opt$external_MAF,".i",opt$internal_MAF,".splmask2.pathway.logistic")
      
    }
    masks <- c("mask1","mask2")
    mask_weight <- c(1, 0.8)
  }else{
    if(opt$variant_type=="syn"){
      if(opt$splmask3 == "splmask3"){
        output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                              ".",subgroup_name, ".", opt$transcript_type,".syn.e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.pathway.logistic")
        
      }else{
        output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                              ".",subgroup_name, ".", opt$transcript_type,".syn.e",opt$external_MAF,".i",opt$internal_MAF,".splmask2.pathway.logistic")
        
      }
      masks <- c("syn")
      mask_weight <- c(1)
    }else{
      if(opt$variant_type=="lof_hc"){
        if(opt$splmask3 == "splmask3"){
          output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                                ".",subgroup_name, ".", opt$transcript_type,".lof_hc.e",opt$external_MAF,".i",opt$internal_MAF,".splmask3.pathway.logistic")
          
        }else{
          output.file <- paste0(dirname(opt$path_name), "/../subgroup/", gsub(".RData", "", basename(opt$pathway_file)), "/", subgroup_name, "/", basename(opt$path_name),
                                ".",subgroup_name, ".", opt$transcript_type,".lof_hc.e",opt$external_MAF,".i",opt$internal_MAF,".splmask2.pathway.logistic")
          
        }
        masks <- c("mask1")
        mask_weight <- c(1)
      }
    }
  }
}

print(opt)
print(input.file)
print(output.file)

## 
sam.info <- read.table(opt$label_file, sep = "\t", stringsAsFactors = F)
colnames(sam.info) <- c("sam_name", "class")
sam.control <- sam.info$sam_name[sam.info$class=="control"]
sam.case <- sam.info$sam_name[sam.info$class=="case"]
print(paste0("the number of subgroup sample: ", length(subgroup_sam)))
print(paste0("the number of case sample: ", length(sam.case)))
sam.case <- intersect(sam.case, subgroup_sam)
print(paste0("the number of subgroup & case sample: ", length(sam.case)))
sam.n.ctrl <- length(sam.control)
sam.n.case <- length(sam.case)
sam.info <- data.frame(sam_name = c(sam.case, sam.control),
                       class = c(rep("case", sam.n.case), rep("control", sam.n.ctrl)),
                       stringsAsFactors = F
)
case.idx <- sam.info$class=="case"
cont.idx <- sam.info$class=="control"

# ，
vrt2sam <- read.table(input.file, sep = "\t", stringsAsFactors = F, skip = 1) 
vrt2sam <- vrt2sam[vrt2sam[,2]%in%masks,] # masks

geneInfo.h <- hash()
for(i in 1:nrow(vrt2sam)){
  ln.ele <- vrt2sam[i,]
  ln.mask <- ln.ele[2]
  ln.gene <- str_split(ln.ele[3], pattern = ";") %>% unlist # 
  ln.mtt.sam <- str_split(ln.ele[4], pattern = ";") %>% unlist # 
  ln.mtt.gnt <- str_split(ln.ele[6], pattern = ";") %>% unlist
  ln.mtt.gnt <- ifelse(grepl("^1/1", ln.mtt.gnt), 2, 1)
  ln.mis.sam <-  str_split(ln.ele[7], pattern = ";") %>% unlist # 
  sapply(ln.gene, function(x){
    if(!has.key(x, geneInfo.h)){
      hash:::.set(geneInfo.h, keys = x,
                  values = list(list(vrt = ln.ele[1], 
                                     vrt.mask = ln.mask,
                                     mtt.sam = ln.mtt.sam,
                                     mtt.gnt = ln.mtt.gnt,
                                     mis.sam = ln.mis.sam
                  )))
    }else{
      hash:::.set(geneInfo.h, keys = x, 
                  values = c(values(geneInfo.h, x), 
                             list(list(vrt = ln.ele[1],
                                       vrt.mask = ln.mask,
                                       mtt.sam = ln.mtt.sam,
                                       mtt.gnt = ln.mtt.gnt,
                                       mis.sam = ln.mis.sam))))
    }
  })
}

# pathway
load(opt$pathway_file)

y.b <- ifelse(sam.info$class=="control", 0, 1)
# 
Z.all <- c()
Z.gene.all <- c()
for(i in keys(geneInfo.h)){ # 
  Z = matrix(0, nrow = length(y.b), ncol = length(geneInfo.h[[i]]))
  rownames(Z) = sam.info$sam_name
  # mask
  vrt.masks <- c()
  for(j in 1:length(geneInfo.h[[i]])){ # ,mask，0，1，2，NAgenotype
    vrt.masks <- c(vrt.masks,
                   geneInfo.h[[i]][[j]]$vrt.mask
    ) %>% unlist() # mask
    Z[geneInfo.h[[i]][[j]]$mtt.sam,j] <- 
      geneInfo.h[[i]][[j]]$mtt.gnt # 0，1，2
  }
  
  vrt.masks.weights <- mask_weight[match(vrt.masks, masks)]
  Z[Z == 2] <- 1 # 012 -> 01
  if(ncol(Z) == 1){ # 
    Z <- apply(Z, 1, function(x){x*vrt.masks.weights}) %>% as.matrix()
  }else{
    Z <- apply(Z, 1, function(x){x*vrt.masks.weights}) %>%
      t() %>% apply(., 1, max)
  }

  Z.all <- cbind(Z.all, Z) # ，
  Z.gene.all <- c(Z.gene.all, i)
}


# pathwaygenotyep
p.result <- c()
# geneInfo.h <- as.list(geneInfo.h)
for(p in 1:length(pathway)){ # 
  path.genes <- pathway[[p]] %>% unique()
  path.genes <- path.genes[path.genes %in% Z.gene.all] # burden
  if(length(path.genes) == 0){next} 
  Z.path <- Z.all[,Z.gene.all %in% path.genes,drop = F] # ，，，012
  Z <- Z.path
  
  path.length <- length(path.genes)
  
  # 
  sam.weighted.masks <- apply(Z, 1, sum)
  
  case.weighted.masks.mean <- sam.weighted.masks[case.idx] %>% mean() %>% round(., digits = 3)
  cont.weighted.masks.mean <- sam.weighted.masks[cont.idx] %>% mean() %>% round(., digits = 3)
  
  case.num <- case.weighted.masks.mean # paste0(case.mut.num,"/",length(case.masks)) # 
  cont.num <- cont.weighted.masks.mean # paste0(cont.mut.num,"/",length(cont.masks)) # 
  
  # # 
  # Z <- apply(Z,2,function(x){
  #   # toFill <- 2*sum(x, na.rm = T)/(sum(!is.na(x))*2) # "fixed" imputes missing genotypes by assigning the mean genotype values (2p).
  #   toFill <- 0 # imputes missing genotypes by 0.
  #   x[is.na(x)] <- toFill
  #   return(x)
  # })
  
  # pathway，
  Z.sum <- sam.weighted.masks
  
  # logistic
  glm.fit <- glm(y.b ~ Z.sum, family = binomial)
  glm.summary <- summary(glm.fit)
  ## 
  beta <- glm.summary$coefficients["Z.sum", "Estimate"]
  std.e <- glm.summary$coefficients["Z.sum", "Std. Error"]
  pValue <- glm.summary$coefficients["Z.sum", "Pr(>|z|)"]
  OR <- exp(beta) %>% signif(., digits = 3)
  OR_CI_L <- exp(beta - 1.96*std.e) %>% signif(., digits = 3)
  OR_CI_U <- exp(beta + 1.96*std.e) %>% signif(., digits = 3)
  
  p.result <- rbind(p.result, c(gene = names(pathway)[p], 
                                path.length,
                                case.num,
                                cont.num,
                                or = OR,
                                ci = paste0(OR_CI_L,"-",OR_CI_U),
                                or_p = pValue,
                                beta = beta,
                                SE = std.e
                                
  ))
}

# p.adjust
p.adj <- apply(p.result[,c(7), drop = F], 2, function(x){
  x <- as.numeric(x)
  x <- p.adjust(x, method = "fdr")
  return(x)
})
p.result <- cbind(p.result, p.adj)

colnames(p.result) <- c("Pathway",
                        "Gene count",
                        "case.gN.mean",
                        "cont.gN.mean",
                        "OR",
                        "95% CI",
                        "OR-P",
                        "BETA",
                        "SE",
                        "P.ADJ"
)

p.result <- p.result[order(as.numeric(p.result[,"OR-P"])),]
write.table(p.result, file = output.file, row.names = F, col.names = T, sep = "\t", quote = F)
date()
