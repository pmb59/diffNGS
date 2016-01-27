# diffNGS
<h4> An R script to identify differential regions in normalized bigwig datasets of ChIP-seq, ATAC-seq, etc. </h4>


<h5>Contact:</h5> 
pm12 [AT] sanger.ac.uk



<h5> Processed GTF: </h5>



    $ cat ensembl_76_transcriptome-GRCh38_15.gtf  | awk '$3 == "gene" {print $0}' > genes_ensembl_76_transcriptome-GRCh38_15.gtf
    
