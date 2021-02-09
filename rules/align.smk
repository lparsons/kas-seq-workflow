rule bwa_mem:
    input:
        reads=get_fq,
        index_amb="{genome}.amb".format(genome=config["ref"]["fasta"]),
        index_ann="{genome}.ann".format(genome=config["ref"]["fasta"]),
        index_bwt="{genome}.bwt".format(genome=config["ref"]["fasta"]),
        index_pac="{genome}.pac".format(genome=config["ref"]["fasta"]),
        index_sa="{genome}.sa".format(genome=config["ref"]["fasta"]),
    output:
        "results/bwa_mem/{sample}-{unit}.bam",
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=config["ref"]["fasta"],
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools/picard.
    resources:
        mem="32G",
    threads: 8
    wrapper:
        "0.67.0/bio/bwa/mem"


rule samtools_index:
    input:
        "results/dedup/{sample}-{unit}.bam",
    output:
        "results/dedup/{sample}-{unit}.bam.bai",
    params:
        "",  # optional params string
    wrapper:
        "0.67.0/bio/samtools/index"


rule mark_duplicates:
    input:
        "results/bwa_mem/{sample}-{unit}.bam",
    output:
        bam="results/dedup/{sample}-{unit}.bam",
        metrics="results/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        "REMOVE_DUPLICATES=true",
    resources:
        mem_mb="10000",
        mem="10G",
    wrapper:
        "0.67.0/bio/picard/markduplicates"


rule bwa_index:
    input:
        "{genome}",
    output:
        "{genome}.amb",
        "{genome}.ann",
        "{genome}.bwt",
        "{genome}.pac",
        "{genome}.sa",
    resources:
        mem="32G",
    params:
        prefix="{genome}",
    log:
        "logs/bwa_index/{genome}.log",
    wrapper:
        "0.67.0/bio/bwa/index"
