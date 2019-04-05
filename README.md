# diffNGS
<h4> R scripts to identify differential open chromatin regions in normalized bigwig datasets of replicated ATAC-seq/DNase-seq data or TF ChIP-seq </h4>

Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS. It requires BEDtools and a GTF file with annotated gene features (e.g., from ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens).

This is a modified version of function narrowpeaksDiff.R of the Bioconductor package <a href="http://bioconductor.org/packages/devel/bioc/html/NarrowPeaks.html">
NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA </a>. An application of the package for Arabidopsis datasets is described in <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> Mateos, Madrigal, et al. (2015) **Genome Biology**: 16:31 </a>. 

Please cite NarrowPeaks if you use diffNGS in your research:

<h5>Citation:</h5> 
Madrigal P, Krajewski P. NarrowPeaks: Shape-based Analysis of Variation in ChIP-Seq using Functional PCA. R package version 1.9.4. 2013. <a href="10.18129/B9.bioc.NarrowPeaks"> DOI: 10.18129/B9.bioc.NarrowPeaks </a>.

<h5>Contact:</h5> 
The code is still under active development, and I welcome biology and computational scientists for all kinds of collaborations and communications. Please feel free to contact Dr. Pedro Madrigal at pmb59(AT)cam.ac.uk if you have any question.



