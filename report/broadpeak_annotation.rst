Broadpeaks from MACS2 annotated by [CHiPseeker](http://bioconductor.org/packages/ChIPseeker) [`annotatePeak`](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html#peak-annotation).

All the peak information contained in peakfile is retained. The position and
strand information of nearest genes are reported. The distance from peak to the
TSS of its nearest gene is also reported. The genomic region of the peak is
reported in annotation column. Since some annotation may overlap, ChIPseeker
adopted the following priority in genomic annotation.

* Promoter
* 5’ UTR
* 3’ UTR
* Exon
* Intron
* Downstream
* Intergenic

Downstream is defined as the downstream of gene end.

Extra columns including SYMBOL, GENENAME, ENSEMBL/ENTREZID will be added.

