### Load the package or install it
  if (!require("fda"))      { install.packages("fda") } 
  if (!require("ICSNP"))    { install.packages("ICSNP") } 
### BioC
  source("http://bioconductor.org/biocLite.R")
  if (!require("genomation"))     { biocLite("genomation")    } 
  if (!require("GenomicRanges"))  { biocLite("GenomicRanges") } 
  if (!require("RColorBrewer"))   { biocLite("RColorBrewer")  } 
  
#set working dir
setwd(".../diffNGS")
  
### Label your experimental conditions
CNDS     <- c(rep("Wt",2), rep("Mt",2) )
### Normalized bigwig files
bws      <- c( "Wt_rep1.norm.bw", "Wt_rep2.norm.bw", "Mt_rep1.norm.bw", "Mt_rep2.norm.bw" ) 
### Aggregated list of peak regions after peak calling against a paired control (in BED format)
peaks    <-  "peaks.bed" 
    
### Parameters of the Analysis
adjP <- 1e-6    # Bonferroni adjusted p-value cuttoff
AbsFc <- 2      # absolute fold change cutoff
#Fl <- 1000      # flanks ( Fl bp upstream and dowstream ) around region center to use in the analysis
Nbasis <- 10    # number of B-spline basis used in the functiona PCA analysis
Bins <- 31      # Number of bins in genomation

### Read BED file
bed <- readGeneric(peaks, keep.all.metadata = FALSE)

#Create data.frame structure to store results
results <- data.frame(region.chr=seqnames(bed), region.start=start(bed)-Fl ,region.end=end(bed)+Fl , condition_C2vsC1=  paste(rev(unique(CNDS)), collapse="_vs_"), avg.C2= rep(NA,length(bed)), avg.C1= rep(NA,length(bed)), pval=rep(NA,length(bed)),  Bonferroni.pval=rep(NA,length(bed)),  fc= rep(NA,length(bed)) ,  diff= rep(FALSE,length(bed)), diff.score= rep(NA,length(bed))     )

# Run diffNGS
source("diffNGS.R")  
x <- diffNGS(bedFile= peaks , headerBed=FALSE, bigwigs=bws, conditions=CNDS, pcs = 2, variation = 0.01, nbasis=Nbasis, NB=Bins)

source("processDiffNGS.R")
processDiffNGS(rawPvals = x$p.values, results=results, CNDS = CNDS, AbsFc = AbsFc, adjP = adjP, plotpdf=TRUE, Xmin=-5, Xmax=5, GTF="genes_ensembl_76_transcriptome-GRCh38_15.gtf", W=5e3 )
#Processed GTF:
#cat ensembl_76_transcriptome-GRCh38_15.gtf  | awk '$3 == "gene" {print $0}' > genes_ensembl_76_transcriptome-GRCh38_15.gtf
    


