setwd("/Users/hqzhao/Desktop/burden/cs_bosheng_1006/res/PCA_plot")
library(ggsci)
my.color <- pal_npg()(7)

pdf("./population_PCA.pdf", height = 3.5,width = 7.5)
par(mar = c(4,4,3,3), cex = 1.0,
    cex.main = 1, cex.axis = 1, cex.lab = 1.2, mfrow = c(1,2), mgp = c(2,0.6,0))
# dat <- read.table("./data/cs_bosheng_1006_clean_pca.eigenvec", sep = " ", header = F, stringsAsFactors = F)
# eigenval <- read.table("./data/cs_bosheng_1006_clean_pca.eigenval")
# variance <- eigenval[,1]/sum(eigenval[,1])*100
# 
# clas <- read.table("./data/sam_class_afFlt.txt", stringsAsFactors = F)
# case.sam <- clas$V1[clas$V2 == "case"]
# label <- ifelse(dat$V1 %in% case.sam, my.color[4], my.color[1])
# plot(dat$V3, dat$V4, 
#      #xlim = c(-0.1,0.1), ylim = c(-0.005,0.01),
#      type = 'n',
#      main = 'A',
#      adj = 0.5,
#      xlab = paste0('PC1', " (",round(variance[1],digits = 1),"%)"),
#      ylab = paste0('PC2', " (",round(variance[2],digits = 1),"%)"),
#      font = 1,
#      font.lab = 1)
# points(dat$V3, dat$V4, col = label, pch = 1, cex = 0.5)
# legend('topright',
#        bty = 'n',
#        cex = 1.0,
#        title = '',
#        c("Case","Control"), #c('African', 'Ad Mixed American', 'East Asian','European', 'South Asian'),
#        fill = my.color[c(1,4)])
# 
# dat <- read.table("./data/cs_bosheng_1006_clean_pca_after_samFlt.eigenvec", sep = " ", header = F, stringsAsFactors = F)
# eigenval <- read.table("./data/cs_bosheng_1006_clean_pca_after_samFlt.eigenval")
# variance <- eigenval[,1]/sum(eigenval[,1])*100
# 
# clas <- read.table("./data/sam_class_afFlt.txt", stringsAsFactors = F)
# case.sam <- clas$V1[clas$V2 == "case"]
# label <- ifelse(dat$V1 %in% case.sam, my.color[4], my.color[1])
# plot(dat$V3, dat$V4, 
#      #xlim = c(-0.1,0.1), ylim = c(-0.005,0.01),
#      type = 'n',
#      main = 'B',
#      adj = 0.5,
#      xlab = paste0('PC1', " (",round(variance[1],digits = 1),"%)"),
#      ylab = paste0('PC2', " (",round(variance[2],digits = 1),"%)"),
#      font = 1,
#      font.lab = 1)
# points(dat$V3, dat$V4, col = label, pch = 1, cex = 0.5)
# legend('topright',
#        bty = 'n',
#        cex = 1.0,
#        title = '',
#        c("Case","Control"), #c('African', 'Ad Mixed American', 'East Asian','European', 'South Asian'),
#        fill = my.color[c(1,4)])

# 1000 G
eigenvec <- read.table("./data/cs_bosheng_1006_1000g_pca_after_samFlt.eigenvec", sep = " ", header = F)
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('Principal Component ', c(1:10), sep = '')
eigenval <- read.table("./data/cs_bosheng_1006_1000g_pca_after_samFlt.eigenval")
variance <- eigenval[,1]/sum(eigenval[,1])*100

# read in the PED data
load("/Users/hqzhao/Desktop/burden/cs_bosheng_1006/saminfo/subgroup.RData")
clas <- read.table("./data/sam_class_afFlt.txt", stringsAsFactors = F)
clas <- rbind(clas[clas$V2 == "control",],
              clas[clas$V1 %in% rownames(subgroup)[subgroup$cvm_undiagnosed==1],])
dat <- read.table("./data/cs_bosheng_1006_pca.eigenvec", sep = " ", header = F, stringsAsFactors = F)
dat <- dat[dat$V1 %in% clas$V1,]
PED <- read.table('./data/20130606_g1k.ped', header = TRUE, skip = 0, sep = '\t')
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[,c("Individual.ID","Population")]
case.sam <- clas$V1[clas$V2 == "case"]
PED <- rbind(PED, data.frame(Individual.ID = dat$V2,
                             Population = ifelse(!dat$V1 %in% case.sam, "Control", "Case")))
eigenvec <- eigenvec[rownames(eigenvec) %in% PED$Individual.ID,]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID==rownames(eigenvec))

# set colours
require('RColorBrewer')

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU",
  "Case",
  "Control"))

# col <- colorRampPalette(c(
#   "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
#   "forestgreen","forestgreen","forestgreen","forestgreen",
#   "grey","grey","grey","grey","grey",
#   "royalblue","royalblue","royalblue","royalblue","royalblue",
#   "black","black","black","black","black",
#   "red",
#   "blue"))(length(unique(PED$Population)))[factor(PED$Population)]

col <- colorRampPalette(c(
  rep(my.color[2],7),
  rep(my.color[3],4),
  rep(my.color[5],5),
  rep(my.color[6],5),
  rep(my.color[7],5),
  my.color[1],
  my.color[4]))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)

plot(project.pca[,1], project.pca[,2],
     type = 'n',
     main = 'A',
     adj = 0.5,
     xlab = paste0('PC1', " (",round(variance[1],digits = 1),"%)"),
     ylab = paste0('PC2', " (",round(variance[2],digits = 1),"%)"),
     font = 1,
     font.lab = 1)
points(project.pca[,1], project.pca[,2], col = col, pch = 1, cex = 0.5)
legend('topright',
       bty = 'n',
       cex = 1.0,
       title = '',
       c("AFR","AMR","EAS","EUR","SAS","Case","Control"), #c('African', 'Ad Mixed American', 'East Asian','European', 'South Asian'),
       fill = my.color[c(2,3,5,6,7,1,4)])


plot(project.pca[,1], project.pca[,3],
     type="n",
     main="B",
     adj=0.5,
     xlab = paste0('PC1', " (",round(variance[1],digits = 1),"%)"),
     ylab = paste0('PC3', " (",round(variance[3],digits = 1),"%)"),
     font=1,
     font.lab = 1)
points(project.pca[,1], project.pca[,3], col=col, pch = 1, cex = 0.5)
legend('topright',
       bty = 'n',
       cex = 1.0,
       title = '',
       c("AFR","AMR","EAS","EUR","SAS","Case","Control"), #c('African', 'Ad Mixed American', 'East Asian','European', 'South Asian'),
       fill = my.color[c(2,3,5,6,7,1,4)])
dev.off()

library(Rtsne)
library("scatterplot3d")
#tsne
set.seed(888)
tsne = Rtsne(project.pca, dims = 3, perplexity = 50,initial_dims = 50)# Run TSNE
scatterplot3d(tsne$Y[,1:3], color = col, xlab = "tSNE 1", ylab = "tSNE 2",
              zlab = "tSNE 3", pch = 16, cex.symbols = 0.3, angle = 1000)
