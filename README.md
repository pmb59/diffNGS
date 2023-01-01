[![GitHub Super-Linter](https://github.com/pmb59/diffNGS/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)

# diffNGS: differential peak analysis in Next-Generation Sequencing data with replicates

diffNGS is a modified version of function `narrowpeaksDiff.R` of the Bioconductor package <a href="http://bioconductor.org/packages/3.10/bioc/html/NarrowPeaks.html">
NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA </a>. 

The repository contains an R script that can be used to identify differential read-enriched regions in **normalized bigwig** files obtained from **replicated** ATAC-seq/DNase-seq/ChIP-seq data. Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS. 

<h2>Examples</h2> 
<h4>for ATAC-seq</h4> 
https://github.com/pmb59/endoderm/tree/master/atacseq/diffNGS
<h4>for ChIP-seq</h4> 
An application of the package for Arabidopsis datasets is described <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> here</a>. 

<h2>Citation</h2> 
Please cite this paper if you use diffNGS in your research:
<br/><br/>

Mateos JL, Madrigal P, et al. <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> Combinatorial activities of SHORT VEGETATIVE PHASE and FLOWERING LOCUS C define distinct modes of flowering regulation in Arabidopsis </a>.  **Genome Biology** 16, 31 (2015)



