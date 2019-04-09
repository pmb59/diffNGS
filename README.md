# diffNGS: differential peak analysis in Next-Generation Sequencing data with replicates

diffNGS is a modified version of function narrowpeaksDiff.R of the Bioconductor package <a href="http://bioconductor.org/packages/devel/bioc/html/NarrowPeaks.html">
NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA </a>. An application of the package for Arabidopsis datasets is described in <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> Mateos, Madrigal, et al. (2015) **Genome Biology**: 16:31 </a>. 

The repository contains an R script that can be used to identify differential read-enriched regions in **normalized bigwig** files obtained from **replicated** ATAC-seq/DNase-seq/ChIP-seq data. Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS. 



<h3>Citation:</h3> 
Please cite NarrowPeaks if you use diffNGS in your research:


Madrigal P, Krajewski P. NarrowPeaks: Shape-based Analysis of Variation in ChIP-Seq using Functional PCA. R package version 1.26.0. 2019. <a href="http://bioconductor.org/packages/NarrowPeaks/"> http://bioconductor.org/packages/NarrowPeaks/ DOI: 10.18129/B9.bioc.NarrowPeaks </a>.

<h3>Contact:</h3> 
Please feel free to contact Dr. Pedro Madrigal at pmb59(AT)cam.ac.uk if you have any question.



