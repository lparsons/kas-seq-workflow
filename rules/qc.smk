def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        if is_single_end(wildcards.sample, wildcards.unit):
            return [ f"{u.fq1}" ]
        else:
            return [ f"{u.fq1}", f"{u.fq2}" ]

    else:
        # yes trimming, use trimmed data
        if not is_single_end(wildcards.sample, wildcards.unit):
            # paired-end sample
            return expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], sample=wildcards.sample, unit=wildcards.unit)
        # single end sample
        return [ "trimmed/{sample}-{unit}.fastq.gz".format(sample=wildcards.sample, unit=wildcards.unit) ]


def get_fastq_group(wildcards):
    fastq_files = get_fq(wildcards)
    return fastq_files[int(wildcards.group)-1]


def get_fastqc(wildcards):
    fastqc_files = list()
    for unit in units.itertuples():
        if is_single_end(unit.sample, unit.unit):
            fastqc_files.append("qc/fastqc/{sample}-{unit}_fastqc.zip".format(sample=unit.sample, unit=unit.unit))
        else:
            fastqc_files.extend(
                ["qc/fastqc/{sample}-{unit}.1_fastqc.zip".format(sample=unit.sample, unit=unit.unit),
                 "qc/fastqc/{sample}-{unit}.2_fastqc.zip".format(sample=unit.sample, unit=unit.unit)])
    return(fastqc_files)


rule fastqc:
    input:
        get_fastq_group
    output:
        html="qc/fastqc/{sample}-{unit}.{group}.html",
        zip="qc/fastqc/{sample}-{unit}.{group}_fastqc.zip"
    params:
        config["params"]["fastqc"]
    threads:
        1
    wrapper:
        "0.67.0/bio/fastqc"


rule samtools_stats:
    input:
        "bwa_mem/{sample}-{unit}.bam"
    output:
        "samtools_stats/{sample}-{unit}.txt"
    params:
        extra="",                       # Optional: extra arguments.
    log:
        "logs/samtools_stats/{sample}-{unit}.log"
    wrapper:
        "0.67.0/bio/samtools/stats"


rule deeptools_plot_coverage:
    input:
        bamfiles=expand("dedup/{unit.sample}-{unit.unit}.bam", unit=units.itertuples()),
        bamindicies=expand("dedup/{unit.sample}-{unit.unit}.bam.bai", unit=units.itertuples()),
    output:
        plotfile="deeptools/plot_coverage.png",
        coverage="deeptools/plot_coverage.tab"
    params:
        ""
    threads:
        8
    log:
        "logs/deeptools/plot_coverage.log"
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        plotCoverage \
                --bamfiles {input.bamfiles:q} \
                --plotFile {output.plotfile:q} \
                --outRawCounts {output.coverage:q} \
                --numberOfProcessors {threads} \
                --smartLabels \
                --plotTitle "Read Coverage (deduplicated)" \
                {params} \
                > {log:q} 2>&1 
        """

rule deeptools_multibamsummary:
    input:
        bamfiles=expand("dedup/{unit.sample}-{unit.unit}.bam", unit=units.itertuples()),
        bamindicies=expand("dedup/{unit.sample}-{unit.unit}.bam.bai", unit=units.itertuples()),
    output:
        multibamsummary="deeptools/multibamsummary.npz",
    params:
        ""
    threads:
        8
    log:
        "logs/deeptools/multibamsummary.log"
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        multiBamSummary \
                bins \
                --bamfiles {input.bamfiles:q} \
                --outFileName {output.multibamsummary:q} \
                --numberOfProcessors {threads} \
                --smartLabels \
                {params} \
                > {log:q} 2>&1
        """


rule deeptools_plot_correlation:
    input:
        cordata="deeptools/multibamsummary.npz",
    output:
        plotfile="deeptools/plot_correlation.png",
        cormatrix="deeptools/plot_correlation_matrix.tsv",
    params:
        ""
    log:
        "logs/deeptools/plot_correlation.log"
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        plotCorrelation \
                --corData {input.cordata:q} \
                --corMethod spearman \
                --whatToPlot heatmap \
                --plotFile {output.plotfile:q} \
                --plotTitle "Alignment correlation (Spearman)" \
                --outFileCorMatrix {output.cormatrix:q} \
                {params} \
                > {log:q} 2>&1
        """

rule deeptools_bampefragmentsize:
    input:
        bamfiles=expand("dedup/{unit.sample}-{unit.unit}.bam", unit=units.itertuples()),
        bamindicies=expand("dedup/{unit.sample}-{unit.unit}.bam.bai", unit=units.itertuples()),
    output:
        histogram="deeptools/bampefragmentsize_histogram.png",
        table="deeptools/bampefragmentsize_table.tsv",
        rawfragmentlengths="deeptools/bampefragmentsize_rawfragmentlengths.tsv",
    params:
        ""
    threads:
        8
    log:
        "logs/deeptools/plot_correlation.log"
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        bamPEFragmentSize \
                --bamfiles {input.bamfiles:q} \
                --histogram {output.histogram:q} \
                --plotTitle "Paired End Fragment Size (deduplicated)" \
                --numberOfProcessors {threads} \
                --table {output.table:q} \
                --outRawFragmentLengths {output.rawfragmentlengths:q} \
                {params} \
                > {log:q} 2>&1
        """


#rule faToTwoBit_fa:
#    input:
#        "{genome}.fa"
#    output:
#        "{genome}.fa.2bit"
#    log:
#        "logs/{genome}.fa_to_2bit.log"
#    params:
#        "" # optional params string
#    wrapper:
#        "0.67.0/bio/ucsc/faToTwoBit"
#
#
#rule deeptools_computegcbias:
#    input:
#        bam="dedup/{sample}-{unit}.bam",
#        bai="dedup/{sample}-{unit}.bam.bai",
#        genome_twobit="{genome}.2bit".format(genome=config["ref"]["fasta"]),
#    output:
#        freq="deeptools/gcbias/{sample}-{unit}_freq.txt",
#        plot="deeptools/gcbias/{sample}-{unit}_plot.png",
#    params:
#        ""
#    threads:
#        8
#    log:
#        "logs/deeptools/gcbias/{sample}-{unit}.log",
#    conda:
#        "../../envs/deeptools.yml"


rule plot_fingerprint:
    input:
        bam_files=expand("dedup/{unit.sample}-{unit.unit}.bam", unit=units.itertuples()),
        bam_idx=expand("dedup/{unit.sample}-{unit.unit}.bam.bai", unit=units.itertuples()),
    output:
        # Please note that --plotFile and --outRawCounts are exclusively defined via output files.
        # Usable output variables, their extensions and which option they implicitly call are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/deeptools/plotfingerprint.html.
        fingerprint="deeptools/plot_fingerprint.png",  # required
        # optional output
        counts="deeptools/plot_fingerprint_raw_counts.tab",
        qc_metrics="deeptools/plot_fingerprint_qc_metrics.txt"
    log:
        "logs/deeptools/plot_fingerprint.log"
    params:
        # optional parameters
        ""
    threads:
        8
    wrapper:
        "0.67.0/bio/deeptools/plotfingerprint"


rule multiqc:
    input:
        get_fastqc,
        expand("trimmed/{unit.sample}-{unit.unit}.qc.txt", unit=units.itertuples()),
        expand("samtools_stats/{unit.sample}-{unit.unit}.txt", unit=units.itertuples()),
        expand("dedup/{unit.sample}-{unit.unit}.metrics.txt", unit=units.itertuples()),
        expand("results/macs2/{unit.sample}-{unit.unit}_peaks.xls", unit=units.itertuples()),
        "deeptools/plot_coverage.tab",
        "logs/deeptools/plot_coverage.log",
        "deeptools/plot_correlation_matrix.tsv",
        "deeptools/bampefragmentsize_table.tsv",
        "deeptools/bampefragmentsize_rawfragmentlengths.tsv",
        "deeptools/plot_fingerprint_qc_metrics.txt",
        "deeptools/plot_fingerprint_raw_counts.tab",
        "results/plot_profile/data.tab",
    output:
        "qc/multiqc_report.html"
    params:
        "--config multiqc_config.yaml"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.67.0/bio/multiqc"