# Combined Human and HBV genome

1. Combined fasta files

    cat ../hg38/genome/fasta/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        ../hbv/NC_003977.fasta > hg38-nc_003977.fasta

2. Combine gtf file

    cat ../hg38/annotation/Homo_sapiens.GRCh38.101.gtf \
        ../hbv/NC_003977.cleaned.gtf > hg38-nc_003977.gtf
