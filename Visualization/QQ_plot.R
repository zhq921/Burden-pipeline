library(magrittr) # for %>%
library(getopt)
library(data.table)
library(QQperm)
library(qqman)

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

files <- list.files(opt$path_name, pattern = "*.perm.qvalue$", recursive = T, full.names = T)
for(i in 1:length(files)){
  # output file
  output.file1 <- paste0(files[i], ".qqplotperm.pdf")
  # output.file2 <- paste0(files[i], ".qqplotuniform.pdf")
  # plot
  dat <- fread(input = files[i], header = T, data.table = F)
  rownames(dat) <- dat$Gene
  dat <- dat[,c("burden.p1.perm", "burden.p1", "burden.p1.qvalue")]
  
  # ## uniform
  # # 3 qqman
  # pdf(output.file2, height = 4.5, width = 4.2, useDingbats = F)
  # my.pvalues <- dat$burden.p1
  # z = qnorm(my.pvalues / 2)
  # lambda = round(median(z^2, na.rm = TRUE) / qchisq(0.5,1), 2)
  # ## lambda
  # if(grepl("lof_mis", files[i])){
  #   tit = paste0("rare LoF + Dmis (lambda = ", lambda, ")")
  # }
  # if(grepl("[.]lof[.]", files[i])){
  #   tit = paste0("rare LoF (lambda = ", lambda, ")")
  # }
  # if(grepl("syn", files[i])){
  #   tit = paste0("rare synonymous (lambda = ", lambda, ")")
  # } 
  # qq(my.pvalues, pch = 19, main = tit,
  #    col = "#045a8d",
  #    xlab = expression(Expected ~ ~-log[10](italic(p))), 
  #    ylab = expression(Observed ~ ~-log[10](italic(p))),
  #    mgp = c(1.7,0.5,0))
  # for(j in 1){
  #   y <- -log10(dat[j,2])
  #   exp.pvalues<-(rank(my.pvalues, ties.method="first")+.5)/(length(my.pvalues)+1)
  #   text(x = -log10(exp.pvalues[1]), y = y, labels = rownames(dat)[j], adj = c(0.5,0.5), cex = 1, offset = 2, font = 3)
  # }
  # dev.off()

  
  # 4 qqperm (recommended)
  pdf(output.file1, height = 4.5, width = 4.2, useDingbats = F)
  
  sig <- (dat$burden.p1.qvalue < 0.05)
  sig.num <- sum(sig)
  if(sig.num >= 5){sig.num <- 5}
  if(sig.num == 0){sig.num <- 1}
  sig.col <- ifelse(sig, "#900b1d", "#045a8d")
  
  e <- -log10(dat[,1])
  o <- -log10(dat[,2])
  adjust.xy = T
  if (adjust.xy) {
    xlim <- c(0, (max(e)))
    ylim <- c(0, (max(o)))
  }else {
    max_lim <- max(max(e), max(o)) + 1
    xlim <- c(0, max_lim)
    ylim <- c(0, max_lim)
  }
  plot(x = e, y = o,pch = 19, xlim = xlim, ylim = ylim, 
       col = sig.col,
       xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
       mgp = c(1.7,0.5,0))
  abline(0, 1, col = "black")
  lambd <- QQperm::estlambda2(dat[,1],dat[,2],plot = F)
  ## lambda
  if(grepl("lof_mis", files[i])){
    title(paste0("rare LoF + Dmis (lambda = ", round(lambd$estimate, digits = 2), ")"))
  }
  if(grepl("[.]lof[.]", files[i])){
    title(paste0("rare LoF (lambda = ", round(lambd$estimate, digits = 2), ")"))
  }
  if(grepl("syn", files[i])){
    title(paste0("rare synonymous (lambda = ", round(lambd$estimate, digits = 2), ")"))
  } 
  for(j in 1:sig.num){
    x <- -log10(dat[j,1])
    y <- -log10(dat[j,2])
    text(x = x, y = y, labels = rownames(dat)[j], adj = c(1.2,0.5), cex = 0.7, offset = 2, font = 3)
  }
  legend(x = 2, y = 1, 
         legend = c("Significant", "Not significant"), 
         col = c("#900b1d", "#045a8d"), pch = 19, cex = 0.8)
  dev.off()
}