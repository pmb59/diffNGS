# diffNGS: differential peak analysis in Next-Generation Sequencing data

diffNGS is a modified version of function narrowpeaksDiff.R of the Bioconductor package <a href="http://bioconductor.org/packages/devel/bioc/html/NarrowPeaks.html">
NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA </a>. An application of the package for Arabidopsis datasets is described in <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> Mateos, Madrigal, et al. (2015) **Genome Biology**: 16:31 </a>. 

The repository contains an R script that can be used to identify differential read-enriched regions in **normalized bigwig** files obtained from **replicated** ATAC-seq/DNase-seq/ChIP-seq data. Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS. 


```R
### Load the package or install it
### CRAN
if (!require("fda"))      { install.packages("fda") } 
if (!require("ICSNP"))    { install.packages("ICSNP") } 
### Bioconductor
source("http://bioconductor.org/biocLite.R")
if (!require("genomation"))     { biocLite("genomation")    } 
if (!require("GenomicRanges"))  { biocLite("GenomicRanges") } 
if (!require("RColorBrewer"))   { biocLite("RColorBrewer")  } 

#set working dir
setwd(".../diffNGS")

### Label your experimental conditions
CNDS     <- c('S1','S1','S2','S2')
### Normalized bigwig files
bws      <- c('S11.bw','S12.bw','S21.bw','S22.bw')
### Aggregated list of peak regions after peak calling against a paired control (in BED format)
peaks    <-  'regions.bed'

### Parameters of the Analysis
Fl <- 0         # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
Nbasis <- 10    # number of B-spline basis used in the functiona PCA analysis
Bins <- 31      # Number of bins in genomation

### Read BED file
bed <- readGeneric(peaks, keep.all.metadata = FALSE)

#Create data.frame structure to store results
results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse="_vs_"), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )
head(results)

# Run diffNGS
source("diffNGS.R")  
x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis=Nbasis, NB=Bins)
head(x)

#Obtain P.values and Fold-changes
for (j in 1:length(x$p.values)  ){

  results$pval[j] <- x$p.values[[j]][1,2]
  
  results$avg.C2[j] <- mean( max(x$fdaprofiles[[3]][j,]) , max(x$fdaprofiles[[4]][j,]) ) 
  results$avg.C1[j] <- mean( max(x$fdaprofiles[[1]][j,]) , max(x$fdaprofiles[[2]][j,]) ) 
  #get the FCs
  if ( results$avg.C2[j] >= results$avg.C1[j]) {  results$log2fc[j]   <-   log2(  ( results$avg.C2[j] + 0.001) / ( results$avg.C1[j] + 0.001  )  )  }    # Increase
  if ( results$avg.C2[j] <  results$avg.C1[j]) {  results$log2fc[j]   <-   log2(  ( results$avg.C1[j] + 0.001) / ( results$avg.C2[j]  + 0.001 )  )  }    # Decrease
  
}

#as.numeric(noquote(unlist(format(.Machine)))[1])
zeroP <- which(results$pval == 0.0)
results$pval[zeroP] <- as.numeric(noquote(unlist(format(.Machine)))[1])

#Correct P-values for multiple hypothesis testing
results$fdr <- p.adjust(results$pval, method = "fdr")
head(results)

# Sort results
results_sorted  <- results[order(results$fdr, decreasing = FALSE),] 
head(results_sorted )

#write results
write.csv(results_sorted , file='results_sorted.csv', row.names = F)


```


<h3>Citation:</h3> 
Please cite NarrowPeaks if you use diffNGS in your research:


Madrigal P, Krajewski P. NarrowPeaks: Shape-based Analysis of Variation in ChIP-Seq using Functional PCA. R package version 1.9.4. 2013. <a href="http://bioconductor.org/packages/NarrowPeaks/"> http://bioconductor.org/packages/NarrowPeaks/ DOI: 10.18129/B9.bioc.NarrowPeaks </a>.

<h3>Contact:</h3> 
The code is still under active development, and I welcome biology and computational scientists for all kinds of collaborations and communications. Please feel free to contact Dr. Pedro Madrigal at pmb59(AT)cam.ac.uk if you have any question.



