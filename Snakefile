import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")



def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        if is_single_end(**wildcards):
            return { 'fq1': f"{u.fq1}" }
        else:
            return { 'fq1': f"{u.fq1}",
                     'fq2': f"{u.fq2}" }

    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return dict(zip(
                ['fq1', 'fq2' ],
                expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)))
        # single end sample
        return { 'fq1': "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) }

##### target rules #####

rule all:
    input:
        #expand(["results/diffexp/{contrast}.diffexp.tsv",
        #        "results/diffexp/{contrast}.ma-plot.svg"],
        #       contrast=config["diffexp"]["contrasts"]),
        #"results/pca.svg",
        "qc/multiqc_report.html"

##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/trim.smk"
#include: "rules/align.smk"
#include: "rules/diffexp.smk"
include: "rules/qc.smk"

