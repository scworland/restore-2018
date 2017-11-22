
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
