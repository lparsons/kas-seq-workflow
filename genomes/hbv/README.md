# HBV Annotation records

1. The fasta sequence and annotations (gff3) were downloaded from NCBI for
   `NC_003977.2`

2. Convert GFF3 to GTF

    gffread -v -F -G -T NC_003977.gff3 -o NC_003977.gtf
    gffread -v -T --force-exons NC_003977.gff3 -o NC_003977.gtf

3. Remove lines without `gene_id`

    grep -v gene_id NC_003977.gtf > NC_003977.cleaned.gtf
