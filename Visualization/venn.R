# Load library
library(ggvenn)
library(readxl)

setwd("/Users/hqzhao/Desktop/burden/cs_bosheng_1006/res/venn_enrichment/")
burden <- read_excel("./burden_p0.05.xlsx", col_names = F)
burden <- burden$...1
gwas <- read.table("./gwas.p.0.05.txt", header = F, stringsAsFactors = F)
gwas <- gwas$V1

dat <- list(
  Burden = burden,
  GWAS = gwas
)

p <- ggvenn(
  dat, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(p, filename = "./venn.pdf", height = 4, width = 4)

# ### reference
# set.seed(20190708)
# genes <- paste("gene",1:1000,sep="")
# x <- list(
#   A = sample(genes,300), 
#   B = sample(genes,525), 
#   C = sample(genes,440),
#   D = sample(genes,350)
# )
# 
# # Change category names
# # Change the fill color
# names(x) <- c("Stage 1","Stage 2","Stage 3", "Stage4")
# ggvenn(
#   x, 
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4
# )
