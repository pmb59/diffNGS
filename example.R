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

### Label your experimental conditions (S1, S2)
CNDS     <- c('S1','S1','S2','S2')
### Normalized bigwig files
bws      <- c('S11.bw','S12.bw','S21.bw','S22.bw')
### Aggregated list of peak regions after peak calling against a paired control (in BED format)
peaks    <-  'regions.bed'

### Parameters of the Analysis
Fl <- 0         # extend flanks ( Fl bp upstream and dowstream ) around regions
Nbasis <- 10    # number of B-spline basis used in the functiona PCA analysis
Bins <- 31      # number of bins in which each region is divided

### Read BED file
bed <- readGeneric(peaks, keep.all.metadata = FALSE)

# Create data.frame to store results
results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse="_vs_"), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  fdr=rep(NA,length(bed)),  log2fc= rep(NA,length(bed))     )
head(results)

# Run diffNGS for 2 Principal Components
source("diffNGS.R")  
x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.3, nbasis = Nbasis, NB = Bins)
head(x)

# Obtain p.values and fold-changes
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

#Correct p-values for Multiple Hypothesis Testing
results$fdr <- p.adjust(results$pval, method = "fdr")

# Sort results by FDR
results_sorted  <- results[order(results$fdr, decreasing = FALSE),] 
head(results_sorted )

# Write results
write.csv(results_sorted , file="diffNGS_results_sorted.csv", row.names = FALSE)
