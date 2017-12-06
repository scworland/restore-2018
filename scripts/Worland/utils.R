
# remove aux sites from dv list
remove_aux <- function(l){
  ncols <- data.frame(ncols=unlist(lapply(l,function(x){ncol(x)})))
  ncols$siteno <- rownames(ncols)
  aux_siteno <- ncols$siteno[which(ncols$ncols !=10)]
  
  for(i in 1:length(aux_siteno)){
    hold <- l[[aux_siteno[i]]]
    hold <- hold[ , -which(names(hold) %in% c("AUX.GAGE_Flow","AUX.GAGE_Flow_cd"))]
    l[[aux_siteno[i]]] <- hold
  }
  
  return(l)
}

# calculate empirical fdcs
sw_efdc <- function(Q){
  df <- data.frame(Q) %>%
    mutate(n = nrow(.),
           m = min_rank(desc(Q)),
           p = m/(n+1)) 
  
  p <- as.numeric(df$p)
}

# generate quantiles over range based of L-moments
sw_lmoms <- function(x, type="pe3") { 
  moms <- lmomco::lmoms(x, no.stop=TRUE)
  if(!are.lmom.valid(moms)) return(NULL) 
  pars <- lmomco::lmom2par(moms, type=type)
  if(is.null(pars)) return(NULL) 
  return(1-plmomco(x, pars))
}

