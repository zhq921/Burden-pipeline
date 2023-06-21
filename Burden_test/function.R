get_max_af <- function(x){
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  return(max(x))
}

get_variant_type <- function(dat){ # VEP：INFOCSQ，
  vrt.csq <- dat[vrt.csq.idx]
  vrt.lofee <- dat[loftee.idx]
  vrt.symbol <- dat[symbol.idx]
  vrt.spliceai <- dat[spliceai.idx] %>% as.numeric()
  vrt.spliceai.symbol <- dat[spliceai.symbol.idx]
  if(
    (grepl(c("stop_gained|frameshift_variant|start_lost|splice_acceptor_variant|splice_donor_variant"), vrt.csq)) &
    (vrt.lofee=="HC")
    ){return("mask1")}
  if((any(vrt.spliceai>0.5) %in% T) & (vrt.symbol == vrt.spliceai.symbol)){
    return("mask2")}
  if((grepl(c("stop_gained|frameshift_variant|start_lost|splice_acceptor_variant|splice_donor_variant"), vrt.csq))){
    return("mask2")}
  if(grepl(c("synonymous_variant"), vrt.csq)){return("syn")}
  if(!grepl(c("inframe_deletion|inframe_insertion|missense_variant|stop_lost"), vrt.csq)){
    return("t-other")}
  # inframe masks
  if(grepl(c("inframe_deletion|inframe_insertion|stop_lost"), vrt.csq)){
    vrt.cadd <- dat[cadd.idx] %>% as.numeric()
    if(sum(vrt.cadd > 20, na.rm = T)==1){return("mask3")}
    if(sum(vrt.cadd > 10, na.rm = T)==1){return("mask4")}
    return("mask5")
  }
  # mis
  vrt.cadd <- dat[cadd.idx] %>% as.numeric()
  vrt.vest <- dat[vest.idx] %>% as.numeric()
  vrt.dann <- dat[dann.idx] %>% as.numeric()
  vrt.eRaw <- dat[eRaw.idx] %>% as.numeric()
  vrt.ePCraw <- dat[ePCraw.idx] %>% as.numeric()
  vrt.fathmm <- dat[fathmm.idx]
  vrt.fathmm.mkl <- dat[fathmm.mkl.idx]
  vrt.provean <- dat[provean.idx]
  vrt.metasvm <- dat[metasvm.idx]
  vrt.metalr <- dat[metalr.idx]
  vrt.mcap <- dat[mcap.idx] %>% as.numeric()
  vrt.primate <- dat[primate.idx]
  vrt.revel <- dat[revel.idx] %>% as.numeric()
  
  mis.sum <- sum(c(vrt.vest > 0.9,
                   vrt.cadd > 10,
                   vrt.dann > 0.9,
                   vrt.eRaw > 0.9,
                   vrt.ePCraw > 0.9,
                   grepl("D", vrt.fathmm), # ，
                   grepl("D", vrt.fathmm.mkl),
                   grepl("D", vrt.provean),
                   grepl("D", vrt.metasvm),
                   grepl("D", vrt.metalr),
                   vrt.mcap > 0.025,
                   grepl("D", vrt.primate),
                   vrt.revel > 0.5),
                 na.rm = T)
  if(mis.sum >= 10){return("mask3")}
  if(mis.sum >= 8){return("mask4")}
  if(mis.sum >= 6){return("mask5")}
  return("mask6")
}

get_variant_type_revel <- function(dat, spliceMask){ # VEP：INFOCSQ，
  vrt.csq <- dat[vrt.csq.idx]
  vrt.lofee <- dat[loftee.idx]
  vrt.symbol <- dat[symbol.idx]
  vrt.spliceai <- dat[spliceai.idx] %>% as.numeric()
  vrt.spliceai.symbol <- dat[spliceai.symbol.idx]
  if(
    (grepl(c("stop_gained|frameshift_variant|start_lost|splice_acceptor_variant|splice_donor_variant"), vrt.csq)) &
    (vrt.lofee=="HC")
  ){return("mask1")}
  if((any(vrt.spliceai>0.5) %in% T) & (vrt.symbol == vrt.spliceai.symbol)){
    return(spliceMask)}
  if((grepl(c("stop_gained|frameshift_variant|start_lost|splice_acceptor_variant|splice_donor_variant"), vrt.csq))){
    return("mask2")}
  if(grepl(c("synonymous_variant"), vrt.csq)){return("syn")}
  if(!grepl(c("inframe_deletion|inframe_insertion|missense_variant|stop_lost"), vrt.csq)){
    return("t-other")}
  # inframe masks
  if(grepl(c("inframe_deletion|inframe_insertion|stop_lost"), vrt.csq)){
    vrt.cadd <- dat[cadd.idx] %>% as.numeric()
    if(sum(vrt.cadd > 20, na.rm = T)==1){return("mask3")}
    if(sum(vrt.cadd > 10, na.rm = T)==1){return("mask4")}
    return("mask5")
  }
  # mis
  # vrt.revel <- dat[revel.idx] %>% as.numeric() # 20220820
  vrt.revel <- dat[revel.idx] %>% strsplit(., split = "&") %>%
    unlist() %>% as.numeric() %>% max(., na.rm = T) # max；-Inf
  if(is.infinite(vrt.revel)){return("mask6n")} # mis without revel score
  if(vrt.revel ==  1 ){return("mask1m")}
  if(vrt.revel >= 0.8){return("mask2m")}
  if(vrt.revel >= 0.6){return("mask3m")}
  if(vrt.revel >= 0.4){return("mask4m")}
  if(vrt.revel >= 0.2){return("mask5m")}
  if(vrt.revel >=  0 ){return("mask6m")}
}
