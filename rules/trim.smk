def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

ruleorder: cutadapt_pe > cutadapt

rule symlink_fastq_pe:
    input:
        get_fastq
    output:
        fastq1="results/symlink_fastq/{sample}-{unit}.1.fastq.gz",
        fastq2="results/symlink_fastq/{sample}-{unit}.2.fastq.gz",
    log:
        "logs/symlink_fastq/{sample}-{unit}.log"
    shell:
        """
        ln -sfr {input[0]:q} {output.fastq1:q}
        ln -sfr {input[1]:q} {output.fastq2:q}
        """

rule symlink_fastq:
    input:
        get_fastq
    output:
        fastq="results/symlink_fastq/{sample}-{unit}.fastq.gz",
    log:
        "logs/symlink_fastq/{sample}-{unit}.log"
    shell:
        """
        ln -s {input[0]:q} {output.fastq:q} >> {log:q} 2>&1
        """

rule cutadapt_pe:
    input:
        ["results/symlink_fastq/{sample}-{unit}.1.fastq.gz",
         "results/symlink_fastq/{sample}-{unit}.2.fastq.gz"]
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters="-a {} -A {} ".format(
                config["trimming"]["read1-adapter"],
                config["trimming"]["read2-adapter"]),
        others=config["params"]["cutadapt-pe"]
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    threads: 24
    wrapper:
        "0.67.0/bio/cutadapt/pe"


rule cutadapt:
    input:
        "results/symlink_fastq/{sample}-{unit}.fastq.gz"
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["trimming"]["read1-adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    threads: 24
    wrapper:
        "0.67.0/bio/cutadapt/se"

