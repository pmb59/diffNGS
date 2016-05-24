# diffNGS
<h4> R scripts to identify differential open chromatin regions in normalized bigwig datasets of replicated ATAC-seq or DNase-seq data </h4>

Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS. It requires BEDtools and a GTF file with annotated gene features (e.g., from <a href="ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/" Ensembl</a>.

This is a modified version of function narrowpeaksDiff.R of the Bioconductor package <a href="http://bioconductor.org/packages/devel/bioc/html/NarrowPeaks.html">
NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA </a>. An application of the package for Arabidopsis datasets is described in <a href="http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1"> Mateos, Madrigal, et al. (2015) **Genome Biology**: 16:31 </a>. Please cite this paper if you use diffNGS in your research.

<h5>Contact:</h5> 
The code is still under active development, and I welcome biology and computational scientists for all kinds of collaborations and communications. Please feel free to contact Dr. Pedro Madrigal at pmb59(AT)cam.ac.uk if you have any question.



