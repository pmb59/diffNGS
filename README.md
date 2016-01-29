# diffNGS
<h4> An R script to identify differential regions in normalized bigwig datasets of replicated ChIP-seq, ATAC-seq, MeRIP-seq, etc. </h4>


<h5>Contact:</h5> 
pm12 [AT] sanger.ac.uk

<h5> Example: </h5>

    $ CNDS     <- c(rep("Wt",2), rep("Mt",2) )
    $ bws      <- c( "Wt_rep1.norm.bw", "Wt_rep2.norm.bw", "Mt_rep1.norm.bw", "Mt_rep2.norm.bw" )   ### Normalized bigwig files
    $ peaks    <-  "peaks.bed" 
    
    $ ### Parameters of the Analysis
    $ adjP <- 1e-6    # Bonferroni adjusted p-value cuttoff
    $ AbsFc <- 4      # absolute fold change cutoff
    $ Fl <- 1000      # flanks ( bp upstream and dowstream ) around region center to use in the analyses
    $ Nbasis <- 10    # number of B-spline basis used in the functiona PCA analysis
 

<h5> Processed GTF: </h5>

    $ cat ensembl_76_transcriptome-GRCh38_15.gtf  | awk '$3 == "gene" {print $0}' > genes_ensembl_76_transcriptome-GRCh38_15.gtf
    
