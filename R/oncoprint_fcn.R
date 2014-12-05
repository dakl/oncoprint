#' function to create an oncoprint plot
oncoprint <- function(somatic, germline){
  genes <- na.omit(unique(c(as.character(somatic$gene), 
                            as.character(germline$gene))))
  #genes
  
  patients <- na.omit(unique(c(as.character(somatic$patient), 
                               as.character(germline$patient))))
  #patients
  
  somatic$gene <- factor(somatic$gene, levels = genes)
  germline$gene <- factor(germline$gene, levels = genes)
  
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
  
  mutmat <- as.data.frame(matrix(0, nc=length(patients), nr=length(genes)))
  colnames(mutmat) <- patients
  rownames(mutmat) <- genes
  mutmat <- mutmat[rev(gene_y_map$y),]
  mutmat <- incrementMatrix(M=mutmat, events = somatic, inc=1)
  mutmat <- incrementMatrix(M=mutmat, events = germline, inc=2)
  # mutmat <- incrementMatrix(M=mutmat, events = cnv_amp, inc=4)
  # mutmat <- incrementMatrix(M=mutmat, events = cnv_gain, inc=8)
  # mutmat <- incrementMatrix(M=mutmat, events = cnv_homloss, inc=16)
  # mutmat <- incrementMatrix(M=mutmat, events = cnv_hetloss, inc=32)
  # # etc...  
  mutmat <- memoSort(mutmat)
  levels(background.m$patient) <- colnames(mutmat)
  
  somatic$gene_y <- germline$gene_y <- NA
  for(gene in genes){
    idx.somatic  <- which( somatic$gene == gene) 
    idx.germline <- which( germline$gene == gene) 
    if(length(idx.somatic) > 0 & length(idx.germline) > 0){
      somatic$gene_y[idx.somatic] <- gene_y_map$y[which(gene_y_map$gene==gene)] + .25
      germline$gene_y[idx.germline] <- gene_y_map$y[which(gene_y_map$gene==gene)] - .25
    } else {
      somatic$gene_y[idx.somatic] <- gene_y_map$y[which(gene_y_map$gene==gene)] 
      germline$gene_y[idx.germline] <- gene_y_map$y[which(gene_y_map$gene==gene)]
    }    
  }
  
  square_w <- .9
  square_h <- .4
  ggplot(background.m, aes(x=patient, y=gene_y)) + geom_tile(fill="gray", colour="white", size=1.1) + 
    scale_y_continuous(breaks=unique(background.m$gene_y), labels=unique(background.m$gene)) + 
    geom_tile(data=somatic, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=square_w, height=square_h, fill="forestgreen") + 
    geom_tile(data=germline, aes(x=patient, y=gene_y), inherit.aes=FALSE, width=square_w, height=square_h, fill="purple", colour=NA) + 
    theme_bw()
}