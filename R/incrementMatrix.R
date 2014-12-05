#' Function to increment gene matrix
#' @param M data frame, colnames=patients, rownames=genes
#' @param events data frame. Needs columns patient and gene
#' @param inc number by which to increment
incrementMatrix <- function(M, events, inc){
  for(k in 1:nrow(events)){
    g <- as.character(events$gene[k])
    s <- as.character(events$patient[k])
    idx.g <- which(rownames(M)==g)
    idx.s <- which(colnames(M)==s)
    if(!is.na(g)){ 
      M[idx.g,idx.s] <- M[idx.g,idx.s]+inc
    }
  }
  M
}

