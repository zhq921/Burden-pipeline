library(ggplot2)
library(reshape2)

setwd("/Users/hqzhao/Desktop/burden/cs_bosheng_1006/res/square_plot/")
df = read.csv("ewce_burden.csv", header = T)
df$P.value <- -log10(df$p.final)
df$time <- factor(df$time, levels = c("5w", "6w", "7w", "8w", "17w"))

ggplot(df)+
  geom_point(aes(CellType,time,size=P.value,fill=sd.from.mean), shape = 22)+
  scale_fill_distiller(name = "Std. Devs",
                       palette = "RdBu", limit=c(-4,4))+
  scale_size_continuous(name = expression(-log[10](italic(p))), 
                        range=c(1,7), 
                        breaks = c(1, 1.5, 2, 2.5,3), 
                        labels = c("1", "1.5", "2", "2.5","3"))+
  theme_classic()+ #coord_flip()+
  theme(axis.text.x=element_text(angle=45, 
                                 hjust=1, colour = "black", 
                                 family='sans'), 
        axis.text.y=element_text(colour = "black", family='sans'),
        legend.title = element_text(family='sans'),
        plot.margin = margin(b = 0, l = 60, r = 0, t = 28)
        )+
  xlab("")+ylab("")
ggsave(filename = "ewce_burden.pdf",
       height = 4.5, width = 8.2)

# gwas
df = read.csv("ewce_gwas.csv", header = T)
df$P.value <- -log10(df$p.final)
df$time <- factor(df$time, levels = c("5w", "6w", "7w", "8w", "17w"))

ggplot(df)+
  geom_point(aes(CellType,time,size=P.value,fill=sd.from.mean), shape = 22)+
  scale_fill_distiller(name = "Std. Devs",
                       palette = "RdBu", limit=c(-4,4))+
  scale_size_continuous(name = expression(-log[10](italic(p))), 
                        range=c(2,7), 
                        breaks = c(1, 2, 3,4), 
                        labels = c("1", "2", "3","4"))+
  theme_classic()+ 
  theme(axis.text.x=element_text(angle=45, 
                                 hjust=1, colour = "black", 
                                 family='sans'), 
        axis.text.y=element_text(colour = "black", family='sans'),
        legend.title = element_text(family='sans'),
        plot.margin = margin(b = 0, l = 60, r = 0, t = 30)
        )+
  xlab("")+ylab("")
ggsave(filename = "ewce_gwas.pdf",
       height = 4.5, width = 8.2)
