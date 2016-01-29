##USE Proper GTF Files! (same as in RNA-seq)
#cat ensembl_76_transcriptome-GRCh38_15.gtf  | awk '$3 == "gene" {print $0}' > genes_ensembl_76_transcriptome-GRCh38_15.gtf 
#rm ensembl_76_transcriptome-GRCh38_15.gtf
##Count number of gene features 
#wc -l genes_ensembl_76_transcriptome-GRCh38_15.gtf


processDiffNGS <- function(rawPvals, results, CNDS, AbsFc, adjP,  plotpdf=F, Xmin=-5, Xmax=5, GTF=F, W=10e3 ) {
  
  
  # Move results to data.frame and get Fold-Change
  for (j in 1:length(rawPvals)  ){
    
    print(paste("Analyzing peak region #", j, "out of", length(rawPvals), sep=" ") )
    
    results$pval[j] <- rawPvals[[j]][1,2]
    
    results$avg.C2[j] <- mean( max(x$fdaprofiles[[3]][j,]) , max(x$fdaprofiles[[4]][j,]) ) 
    results$avg.C1[j] <- mean( max(x$fdaprofiles[[1]][j,]) , max(x$fdaprofiles[[2]][j,]) ) 
    #get the FCs
    if ( results$avg.C2[j] >= results$avg.C1[j]) {  results$fc[j]   <-    results$avg.C2[j]  / results$avg.C1[j] } # Increase
    if ( results$avg.C2[j] <  results$avg.C1[j]) {  results$fc[j]   <-  - results$avg.C1[j]  / results$avg.C2[j] } # Decrease
    
  }
  
  #as.numeric(noquote(unlist(format(.Machine)))[1])
  zeroP <- which(results$pval == 0.0)
  results$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])
  
  #Correct P-values for multiple hypothesis testing- BH or Bonferroni
  results$Bonferroni.pval <- p.adjust(results$pval, method = "bonferroni")
  #head(results)
  
  #report score ( differential peaks will be ranked by this score )
  results$diff.score <- abs (results$fc ) * -log10(results$Bonferroni.pval)
  
  Down <- which ( results$fc <= -AbsFc  & results$Bonferroni.pval<= adjP )
  if (length(Down) >0 ) { results[Down ,]$diff <- "TRUE.Down"  }
  Up <- which(results$fc >=   AbsFc  & results$Bonferroni.pval<= adjP )
  if (length(Up) >0 ) { results[Up ,]$diff <- "TRUE.Up"  }
  
  #head(results)
  colnames(results) <- c("region.chr",  "region.start",	"region.end",	"condition_C2vsC1",	paste("avg.",rev(unique(CNDS))[1],sep=""),	paste("avg.",rev(unique(CNDS))[2],sep=""),	"pval",	"Bonferroni.pval",	"fc",	"diff", "diff.score" )
  
  #Report all tests in the same order as in the Bed file
  write.csv(results, file=paste(  paste(rev(unique(CNDS)), collapse="_vs_"), "diffNGS.csv",sep="_"), row.names=F)
  
  #Generate density Figure
  
  if(plotpdf == TRUE){
    pdf(file=paste(  paste(rev(unique(CNDS)), collapse="_vs_"), "diffNGS.pdf",sep="_"), height=5, width=5)
    i1 <- which(results$Bonferroni.pval <= adjP)
    i2 <- which(results$Bonferroni.pval > adjP)
    i3 <- which(results$Bonferroni.pval <= adjP & abs(results$fc) >=  AbsFc )
    # Kernel Density Plot
    d1 <- density(results$fc[i1]) # returns the density data
    d2 <- density(results$fc[i2]) # returns the density data
    d3 <- density(results$fc[i3]) # returns the density data
    YL <- max(c(d1$y,d2$y,d3$y))
    plot(frame=F, d1, ylim=c(0,YL), xlab="Peak Mean Fold-change", main=paste(rev(unique(CNDS)), collapse="_vs_") , xlim=c(Xmin,Xmax)) # plots the results 
    points(d2, lty=2, type='l')
    points(d3, lty=1, type='l', lwd=2, col='dodgerblue3')
    abline(v = 1, col="red", lwd=2)
    mtext(paste("total regions = " , nrow(results),"; ", length(Up), " Up; ",length(Down)," Down", sep="") )
    legend("topright", legend=c(expression( P <= adjP ), "n.s.", "diffNGS") , col=c("black","black","dodgerblue3"), bty="n", lty=c(1,2,1), lwd=2)
    dev.off()
    
  }
  
  
  
  results <- results[which(results$diff != FALSE)  ,]
  results <- results[order(results$diff.score, decreasing=TRUE),]
  #head(results)
  #tail(results)
  
  # Include closest genes in the table
  #Target genes (Wkb)
  if ( GTF != F ) {
    ranDomFileName <- paste(sample( LETTERS,20), collapse="")
    write.table(x=results[,1:3], file = ranDomFileName, append = FALSE, quote = F, sep = "\t", row.names = F, col.names = F )
    
    system("chmod + bedtools")
    system( paste("./bedtools window -w", W, "-a",ranDomFileName,"-b", GTF, ">", paste(ranDomFileName,"genes",sep="."), sep=" ") )   #> temp.genes
  
    raw <- read.table(paste(ranDomFileName,"genes",sep="."), head=F)
    
    gene_id   <- c()
    gene_type <- c()
    gene_name <- c()
    
    for (i in 1:nrow(raw)) {
      gene_id[i] <- as.character( raw$V13[i] )  
      gene_type[i] <- as.character(raw$V22[i] )
      gene_name[i] <-as.character(raw$V16[i] ) 
    }
    output <- cbind(raw[,1:3],gene_id,gene_type, gene_name )
    #add gene names, gene ids, and gene_types to 'results'
    gene_ids <- c()
    gene_types <- c()
    gene_names <- c()
    for ( l in 1:nrow(results) ){
      
      l_matched <- which( as.character(output[,1]) ==  results$region.chr[l]  & as.character( output[,2] ) ==  results$region.start[l]  &  output[,3] ==  results$region.end[l]  ) 
      gene_ids[l]   <- paste(output$gene_id[l_matched], collapse=";")
      gene_types[l] <- paste(output$gene_type[l_matched], collapse=";")
      gene_names[l] <- paste(output$gene_name[l_matched], collapse=";")
      
    }
    results <- cbind( results, gene_ids, gene_types, gene_names)
    
    system(paste('rm', ranDomFileName ,sep = " ")   )
    system(paste('rm', paste(ranDomFileName,"genes",sep=".") ,sep = " ")   )
    
  }
  
  
  write.csv(results, file=paste( paste(rev(unique(CNDS)), collapse="_vs_"), "fc",AbsFc, "adjP",adjP, "diffNGS.csv",sep="_"), row.names=F)
  
  write.table(x=results, file = paste( paste(rev(unique(CNDS)), collapse="_vs_"), "fc",AbsFc, "adjP",adjP, "diffNGS.bed",sep="_"), append = FALSE, quote = F, sep = "\t", row.names = F,  col.names = F)
  
}


