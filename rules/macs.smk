rule callpeak_options:
    input:
        treatment="results/dedup/{sample}-{unit}.bam",
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        ### required
        report(
            "results/macs2/{sample}-{unit}_peaks.xls",
            caption="../report/macs2_peaks.rst",
            category="KAS-Seq Peaks",
            subcategory="{sample}-{unit}",
        ),
        ### optional output files
        # these output extensions internally set the --bdg or -B option:
        "results/macs2/{sample}-{unit}_treat_pileup.bdg",
        #"results/macs2_{sample}-{unit}_control_lambda.bdg",
        # these output extensions internally set the --broad option:
        report(
            "results/macs2/{sample}-{unit}_peaks.broadPeak",
            caption="../report/macs2_broadpeak.rst",
            category="KAS-Seq Peaks",
            subcategory="{sample}-{unit}",
        ),
        report(
            "results/macs2/{sample}-{unit}_peaks.gappedPeak",
            caption="../report/macs2_gappedpeak.rst",
            category="KAS-Seq Peaks",
            subcategory="{sample}-{unit}",
        ),
    log:
        "logs/macs2/{sample}-{unit}_callpeak.log",
    params:
        "--format BAM --gsize hs --broad-cutoff 0.01 --qvalue 0.01",
    wrapper:
        "0.67.0/bio/macs2/callpeak"


rule broadpeak_headers:
    input:
        "results/macs2/{sample}-{unit}_peaks.broadPeak",
    output:
        "results/macs2/{sample}-{unit}_peaks.broadPeak.withHeaders",
    log:
        "logs/macs2/{sample}-{unit}_broadpeak_headers.log",
    conda:
        "../envs/coreutils.yml"
    shell:
        """
        echo -e 'chr\tstart\tend\tname\tscore\tstrand\tfold_enrichment\t-log10(pvalue)\t-log10(qvalue)' > {output:q}
        cat {input:q} >> {output:q}
        """
