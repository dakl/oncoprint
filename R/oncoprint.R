#' function to create an oncoprint plot
#'
#' @param M mutation matrix
#' @param keys list with the following elements: somatic, germline, amp, del, upreg, downreg
#' @param sortGenes boolean whether or not to sort the genes, default FALSE
oncoprint <- function(M, keys=list(somatic="MUT", germline="GERMLINE", amp="AMP", 
                      del="HOMDEL", upreg="UP", downreg="DOWN"), sortGenes=FALSE){
  
  #M
  all <- melt(M, varnames = c("gene", "patient"), value.name = "alteration")
  #head(all,n = 50)
  genes <- na.omit(unique(as.character(all$gene)))
  patients <- na.omit(unique( as.character(all$patient) ))
  
  e <- data.frame(gene=NA, patient=NA, alteration=NA)
  data <- list()
  for( alteration in names(keys) ){
    data[[alteration]] <- all[ grep(pattern = keys[[alteration]], all$alteration) ,]
    if(nrow(data[[alteration]]) == 0 ) data[[alteration]] <- e
    data[[alteration]]$gene <- factor(data[[alteration]]$gene, levels=genes)
    data[[alteration]]$gene_y <- 0
  }    
  background <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(background) <- patients
  background$gene <- genes
  
  background.m <- melt(background, id.vars = "gene", variable.name =  "patient")
  
  # http://stackoverflow.com/a/16486873/179444
  # > when you add a new data set to a geom you need to use the 
  # > data= argument. Or put the arguments in the proper order mapping=..., data=.... Take a look at the arguments for ?geom_line.
  #
  # sort genes, by name  
  # 'y' in the df below is the actual y-values used in the plot so it sets the order
  gene_y_map <- data.frame(gene=genes, y=order(genes, decreasing = TRUE))
  background.m$gene_y <- gene_y_map$y[ match( background.m$gene, gene_y_map$gene ) ]
  
  ## create numerical mutation matrix
  mutmat <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(mutmat) <- patients
  rownames(mutmat) <- genes
  mutmat <- mutmat[rev(gene_y_map$y),]
  # from https://github.com/gideonite/WIP/blob/gh-pages/oncoprint/MemoSort.js
  # // sorting order : amplification, deletion, mutation, mrna, rppa
  # // mutation > 0
  # // amp > del > 0
  mutmat <- incrementMatrix(M=mutmat, events = data$amp, inc=128)
  mutmat <- incrementMatrix(M=mutmat, events = data$del, inc=64)
  mutmat <- incrementMatrix(M=mutmat, events = data$somatic, inc=32)
  mutmat <- incrementMatrix(M=mutmat, events = data$germline, inc=16)
  mutmat <- incrementMatrix(M=mutmat, events = data$upreg, inc=8)
  mutmat <- incrementMatrix(M=mutmat, events = data$downreg, inc=4)
  
  mutmat <- memoSort(mutmat, sortGenes = sortGenes)
  levels(background.m$patient) <- colnames(mutmat)
  
  for(gene in genes){
    idx.somatic  <- which( data$somatic$gene == gene) 
    idx.germline <- which( data$germline$gene == gene) 
    if(length(idx.somatic) > 0 & length(idx.germline) > 0){
      data$somatic$gene_y[idx.somatic] <- gene_y_map$y[which(gene_y_map$gene==gene)] + .25
      data$germline$gene_y[idx.germline] <- gene_y_map$y[which(gene_y_map$gene==gene)] - .25
    } else if(length(idx.somatic) > 0) {
      data$somatic$gene_y[idx.somatic] <- gene_y_map$y[which(gene_y_map$gene==gene)] 
    }else if(length(idx.germline) > 0) {
      data$germline$gene_y[idx.germline] <- gene_y_map$y[which(gene_y_map$gene==gene)]
    }        
  }
  data$amp$gene_y <- gene_y_map$y[match(data$amp$gene, gene_y_map$gene)]
  data$del$gene_y <- gene_y_map$y[match(data$del$gene, gene_y_map$gene)]
  data$upreg$gene_y <- gene_y_map$y[match(data$upreg$gene, gene_y_map$gene)]
  data$downreg$gene_y <- gene_y_map$y[match(data$downreg$gene, gene_y_map$gene)]
  
  square_w <- .9
  square_h <- .4
  ggplot(background.m, aes(x=patient, y=gene_y)) + geom_tile(fill="gray", colour="white", size=1.1) + 
    scale_y_continuous(breaks=unique(background.m$gene_y), labels=unique(background.m$gene)) + 
    geom_tile(data=data$amp, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=.9, height=.9, fill="firebrick", colour=NA, size=2) + 
    geom_tile(data=data$del, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=.9, height=.9, fill="blue", colour=NA, size=2) + 
    geom_tile(data=data$somatic, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=square_w, height=square_h, fill="forestgreen") + 
    geom_tile(data=data$germline, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=square_w, height=square_h, fill="purple", colour=NA) + 
    geom_tile(data=data$upreg, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=.9, height=.9, fill=NA, colour="firebrick", size=2) + 
    geom_tile(data=data$downreg, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=.9, height=.9, fill=NA, colour="dodgerblue", size=2) + 
    theme_minimal() + xlab("Sample") + ylab("Gene")
}
