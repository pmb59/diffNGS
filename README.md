# diffNGS: differential peak analysis in Next-Generation Sequencing data with replicates

diffNGS is a modified version of function `narrowpeaksDiff.R` of the Bioconductor package [NarrowPeaks: Shape-based Analysis of Variation in ChIP-seq using Functional PCA](http://bioconductor.org/packages/3.10/bioc/html/NarrowPeaks.html).

The repository contains an R script that can be used to identify differential read-enriched regions in **normalized bigwig** files obtained from **replicated** ATAC-seq/DNase-seq/ChIP-seq data. Peak calling at each condition and generation of normalized bigwig tracks of the samples should be done before using diffNGS.

## Examples

#### for ATAC-seq
[https://github.com/pmb59/endoderm/tree/master/atacseq/diffNGS](https://github.com/pmb59/endoderm/tree/master/atacseq/diffNGS)

#### for ChIP-seq
An application of the package for Arabidopsis datasets is described [here](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1).

## Citation
Please cite this paper if you use diffNGS in your research:


Mateos JL, Madrigal P, et al. [Combinatorial activities of SHORT VEGETATIVE PHASE and FLOWERING LOCUS C define distinct modes of flowering regulation in Arabidopsis. **Genome Biology** 16, 31 (2015)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0597-1).
