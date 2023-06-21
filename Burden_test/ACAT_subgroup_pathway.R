
library(SKAT)
library(ACAT)
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
system(command = paste0("mkdir -p ", dirname(opt$path_name), "/../subgroup/ACAT/", subgroup_name)) # 

opt$external_MAF <- as.numeric(opt$external_MAF)
opt$internal_MAF <- as.numeric(opt$internal_MAF)
input.file <- paste0(dirname(opt$path_name), "/../subgroup/qualified_variant/", subgroup_name, "/", basename(opt$path_name),
                     ".",subgroup_name,".",opt$transcript_type,".e",opt$external_MAF,".i",opt$internal_MAF,".txt")
if(opt$variant_type=="lof_mis"){ # lof_mis，lof+misburden
  output.file <- paste0(dirname(opt$path_name), "/../subgroup/ACAT/", subgroup_name, "/", basename(opt$path_name),
                        ".",subgroup_name,".", opt$transcript_type,".lof_mis.e",opt$external_MAF,".i",opt$internal_MAF,".pathway.ACAT")
  masks <- c("mask1","mask2","mask3","mask4","mask5",
             "mask1m","mask2m","mask3m","mask4m","mask5m")
  mask_weight <- c(1, 0.8, 0.5, 0.3, 0.1, 
                   1, 0.8, 0.5, 0.3, 0.1)
}else{
  if(opt$variant_type=="lof"){
    output.file <- paste0(dirname(opt$path_name), "/../subgroup/ACAT/", subgroup_name, "/", basename(opt$path_name),
                          ".",subgroup_name, ".", opt$transcript_type,".lof.e",opt$external_MAF,".i",opt$internal_MAF,".pathway.ACAT")
    masks <- c("mask1","mask2")
    mask_weight <- c(1, 0.8)
  }else{
    if(opt$variant_type=="syn"){
      output.file <- paste0(dirname(opt$path_name), "/../subgroup/ACAT/", subgroup_name, "/", basename(opt$path_name),
                            ".",subgroup_name, ".", opt$transcript_type,".syn.e",opt$external_MAF,".i",opt$internal_MAF,".pathway.ACAT")
      masks <- c("syn")
      mask_weight <- c(1)
    }else{
      if(opt$variant_type=="lof_hc"){
        output.file <- paste0(dirname(opt$path_name), "/../subgroup/ACAT/", subgroup_name, "/", basename(opt$path_name),
                              ".",subgroup_name, ".", opt$transcript_type,".lof_hc.e",opt$external_MAF,".i",opt$internal_MAF,".pathway.ACAT")
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


if(!is.null(opt$covariate_PCA)){
  covMat <- read.table(opt$covariate_PCA, header = T, stringsAsFactors = F, sep = "\t")
  rownames(covMat) <- covMat[,1]
  covMat <- covMat[sam.info$sam_name,-1]
  X <- as.matrix(covMat)
} # cov

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
      hash:::.set(geneInfo.h, keys = x, values = list(list(vrt = ln.ele[1], 
                                                           vrt.mask = ln.mask,
                                                           mtt.sam = ln.mtt.sam,
                                                           mtt.gnt = ln.mtt.gnt,
                                                           mis.sam = ln.mis.sam
      )))
    }else{
      hash:::.set(geneInfo.h, keys = x, values = c(values(geneInfo.h, x), list(list(vrt = ln.ele[1], 
                                                                                    vrt.mask = ln.mask,
                                                                                    mtt.sam = ln.mtt.sam,
                                                                                    mtt.gnt = ln.mtt.gnt,
                                                                                    mis.sam = ln.mis.sam))))
      
    }
  })
}

# pathway
load(opt$pathway_file)


# ACATp
p.result <- c()
y.b <- ifelse(sam.info$class=="control", 0, 1)
if(exists(("X"))){
  obj <- SKAT_Null_Model(y.b ~ X, out_type="D") # PCA，PCA，1000G，
  obj.acat <- NULL_Model(y.b, Z = X)
}else{
  obj <- SKAT_Null_Model(y.b ~ 1, out_type="D")
  obj.acat <- NULL_Model(y.b, Z = NULL)
}
# geneInfo.h <- as.list(geneInfo.h)
for(p in 1:length(pathway)){ # 
  path.genes <- pathway[[p]]
  path.genes <- path.genes[path.genes %in% keys(geneInfo.h)] # burden
  if(length(path.genes) == 0){next} 
  Z.path <- c()
  vrt.masks.path <- c()
  vrt.weights.path <- c()
  for(i in path.genes){ # ，genotype
    # genotype,，
    Z = matrix(0, nrow = length(y.b), ncol = length(geneInfo.h[[i]]))
    rownames(Z) = sam.info$sam_name
    # mask
    vrt.masks <- c()
    for(j in 1:length(geneInfo.h[[i]])){ # ,mask，0，1，2，NAgenotype
      vrt.masks <- c(vrt.masks,
                     geneInfo.h[[i]][[j]]$vrt.mask
      ) %>% unlist() # mask
      Z[geneInfo.h[[i]][[j]]$mtt.sam,j] <- geneInfo.h[[i]][[j]]$mtt.gnt # 0，1，2
      if(!("" %in% geneInfo.h[[i]][[j]]$mis.sam)){ # missing，NA
        Z[geneInfo.h[[i]][[j]]$mis.sam,j] <- NA
      }
    } 
    vrt.weights <- mask_weight[match(vrt.masks, masks)] # 
    Z.path <- cbind(Z.path, Z)
    vrt.masks.path <- c(vrt.masks.path, vrt.masks)
    vrt.weights.path <- c(vrt.weights.path, vrt.weights)
  }
  Z <- Z.path
  vrt.masks <- vrt.masks.path
  vrt.weights <- vrt.weights.path
  mask.out <- paste(table(vrt.masks),names(table(vrt.masks)), sep = "_", collapse = "|")
  # mask
  sam.masks <- apply(Z,1,function(x) {
    tmp <- (x %in% c(1,2)) # idx
    if(!any(tmp)){ # ，
      return("no_mutation")
    }
    vrt.masks[tmp] %>% min() %>% return() # ，mask
  })
  case.masks <- sam.masks[case.idx]
  case.mut.num <- sum(case.masks!="no_mutation")
  case.num <- paste0(case.mut.num,"/",length(case.masks))
  case.masks.out <- grep("no_mutation", case.masks, invert = T, value = T)
  case.masks.out <- paste(table(case.masks.out),names(table(case.masks.out)), sep = "_", collapse = "|") %>% 
    paste0(case.mut.num,"/",length(case.masks),":",.)
  cont.masks <- sam.masks[cont.idx]
  cont.mut.num <- sum(cont.masks!="no_mutation")
  cont.num <- paste0(cont.mut.num,"/",length(cont.masks))
  cont.masks.out <- grep("no_mutation", cont.masks, invert = T, value = T)
  cont.masks.out <- paste(table(cont.masks.out),names(table(cont.masks.out)), sep = "_", collapse = "|") %>% 
    paste0(cont.mut.num,"/",length(cont.masks),":",.)
  fisher.p <- fisher.test(
    matrix(c(case.mut.num, length(case.masks) - case.mut.num,
             cont.mut.num, length(cont.masks) - cont.mut.num), nrow = 2, byrow = T),
    alternative = "greater"
  )$p.value
  
  # 
  Z <- apply(Z,2,function(x){
    # toFill <- 2*sum(x, na.rm = T)/(sum(!is.na(x))*2) # "fixed" imputes missing genotypes by assigning the mean genotype values (2p).
    toFill <- 0 # imputes missing genotypes by 0.
    x[is.na(x)] <- toFill
    return(x)
  })
  
  # OR_calculator
  OR_calculator <- function(dat){
    dat[1,1] <- as.numeric(dat[1,1]) # double，
    or <- (dat[2,2]*dat[1,1])/(dat[1,2]*dat[2,1])
    or.se <- sqrt(1/dat[1] + 1/dat[2] + 1/dat[3] +1/dat[4])
    or.low <- exp(log(or) - 1.96*or.se)
    or.up <- exp(log(or) + 1.96*or.se)
    chisq.stat <- (sum(dat) - 1) * (dat[4]*dat[1] - dat[2]*dat[3])^2 /
      ((dat[1]+dat[3]) * (dat[2]+dat[4]) * (dat[1]+dat[2]) * (dat[3]+dat[4]))
    p <- pchisq(chisq.stat, df = 1, lower.tail = F)
    res <- c(signif(or, digits = 3), 
             paste0(signif(or.low,digits = 3), "-",
                    signif(or.up,digits = 3)),
             signif(p,digits = 3)
    )
    return(res)
  }
  # OR
  get_or <- function(case.masks, cont.masks, masks, mask_weight){ 
    case.weighted.num <- mask_weight[match(case.masks, masks)] %>% sum() # 
    cont.weighted.num <- mask_weight[match(cont.masks, masks)] %>% sum()
    case.weighted.num.nomut <- length(case.masks) - case.weighted.num
    cont.weighted.num.nomut <- length(cont.masks) - cont.weighted.num
    dat <- matrix(c(
      (case.weighted.num+0.5), (case.weighted.num.nomut+0.5),
      (cont.weighted.num+0.5), (cont.weighted.num.nomut+0.5)
    ), byrow = T, nrow = 2)
    res <- OR_calculator(dat = dat)
    # Since we are interested in rare variants and the corresponding 2 × 2 table may consist 
    # of entries with 0 observations, we applied the amended estimator of the odds ratio by 
    # adding 0.5 to each cell. It has been suggested that the amended estimator of the odds ratio behaves well
    return(res)
  }
  or <- get_or(case.masks = case.masks, cont.masks = cont.masks, 
               masks = c(masks,"no_mutation"), mask_weight = c(mask_weight, 0))
  
  # P
  burden.p <- SKAT(Z, obj,r.corr = 1, weights = vrt.weights)$p.value
  skat.p <- SKAT(Z, obj,r.corr = 0, weights = vrt.weights)$p.value
  acatv.p <- ACAT_V(Z, obj.acat, weights = vrt.weights)
  p.result <- rbind(p.result, c(gene = names(pathway)[p], 
                                case.num,
                                cont.num,
                                or = or[1],
                                ci = or[2],
                                or_p = or[3],
                                mask.out,
                                case.masks.out,
                                cont.masks.out,
                                fisher.p,
                                burden.p1 = burden.p[1],
                                # burden.p2 = burden.p[2],
                                skat.p1 = skat.p[1],
                                # skat.p2 = skat.p[2],
                                acatv.p1 = acatv.p[1]# ,
                                # acatv.p2 = acatv.p[2],
                                # skato.p1 = skato.p[1],
                                # skato.p2 = skato.p[2],
                                # acato.p = acato.p
  ))
}

# p.adjust
p.adj <- apply(p.result[,c(10:13)], 2, function(x){
  x <- as.numeric(x)
  x <- p.adjust(x, method = "fdr")
  return(x)
})
p.result <- cbind(p.result, p.adj)

colnames(p.result) <- c("Pathway",
                        "case.num",
                        "cont.num",
                        "OR",
                        "95% CI",
                        "OR-P",
                        "vrt.masks",
                        "case.masks",
                        "cont.masks",
                        "fisher.p",
                        "burden.p1",# "burden.p2", 
                        "skat.p1",# "skat.p2",
                        "acatv.p1",# "acatv.p2",
                        # "skato.p1","skato.p2",
                        # "acato.p",
                        "fisher.p.adj",
                        "burden.p1.adj",# "burden.p2.adj", 
                        "skat.p1.adj",# "skat.p2.adj",
                        "acatv.p1.adj"# , "acatv.p2.adj",
                        # "skato.p1.adj","skato.p2.adj",
                        # "acato.p.adj"
                        )
p.result <- p.result[order(as.numeric(p.result[,"acato.p"])),]
write.table(p.result, file = output.file, row.names = F, col.names = T, sep = "\t", quote = F)
date()

