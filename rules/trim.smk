ruleorder: cutadapt_pe > cutadapt


rule symlink_fastq_pe:
    input:
        unpack(get_fastq),
    output:
        fastq1="results/symlink_fastq/{sample}-{unit}.1.fastq.gz",
        fastq2="results/symlink_fastq/{sample}-{unit}.2.fastq.gz",
    log:
        "logs/symlink_fastq/{sample}-{unit}.log",
    conda:
        "../envs/coreutils.yml"
    shell:
        """
        ln -sfr {input.fq1:q} {output.fastq1:q}
        ln -sfr {input.fq2:q} {output.fastq2:q}
        """


rule symlink_fastq:
    input:
        unpack(get_fastq),
    output:
        fastq="results/symlink_fastq/{sample}-{unit}.fastq.gz",
    log:
        "logs/symlink_fastq/{sample}-{unit}.log",
    conda:
        "../envs/coreutils.yml"
    shell:
        """
        ln -s {input.fq1:q} {output.fastq:q} >> {log:q} 2>&1
        """


rule cutadapt_pe:
    input:
        [
            "results/symlink_fastq/{sample}-{unit}.1.fastq.gz",
            "results/symlink_fastq/{sample}-{unit}.2.fastq.gz",
        ],
    output:
        fastq1="results/trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}-{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        adapters="-a {} -A {} ".format(
            config["trimming"]["read1-adapter"], config["trimming"]["read2-adapter"]
        ),
        others=config["params"]["cutadapt-pe"],
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    threads: 24
    wrapper:
        "0.67.0/bio/cutadapt/pe"


rule cutadapt:
    input:
        "results/symlink_fastq/{sample}-{unit}.fastq.gz",
    output:
        fastq="results/trimmed/{sample}-{unit}.fastq.gz",
        qc="results/trimmed/{sample}-{unit}.qc.txt",
    params:
        "-a {} {}".format(
            config["trimming"]["read1-adapter"], config["params"]["cutadapt-se"]
        ),
    log:
        "logs/cutadapt/{sample}-{unit}.log",
    threads: 24
    wrapper:
        "0.67.0/bio/cutadapt/se"
