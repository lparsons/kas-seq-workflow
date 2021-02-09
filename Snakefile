import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####


configfile: "config.yaml"


validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

##### target rules #####


rule all:
    input:
        #expand(["results/diffexp/{contrast}.diffexp.tsv",
        #        "results/diffexp/{contrast}.ma-plot.svg"],
        #       contrast=config["diffexp"]["contrasts"]),
        #"results/pca.svg",
        expand("results/bwa_mem/{unit.sample}-{unit.unit}.bam", unit=units.itertuples()),
        expand(
            "results/macs2/{unit.sample}-{unit.unit}_peaks.xls", unit=units.itertuples()
        ),
        "results/plot_heatmap/heatmap.png",
        "results/plot_profile/profile.png",
        "results/peak_analysis/broadpeakAnnotationDistributionPlot.pdf",
        "results/peak_analysis/broadPeakAnnoList.Rdata",
        expand(
            "results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.tsv.gz",
            unit=units.itertuples(),
        ),
        expand(
            "results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.summary.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/peak_analysis/{unit.sample}-{unit.unit}_broadpeak_overlap_enrichment.tsv",
            unit=units.itertuples(),
        ),
        "results/peak_analysis/broadPeakAvgProfilePlot.pdf",
        "results/peak_analysis/broadPeakHeatmapPlot.pdf",
        "results/peak_analysis/broadPeakMatrixList.Rdata",
        "results/peak_analysis/broadPeakAnnotationVennDiagram.pdf",
        "qc/multiqc_report.html",


##### setup singularity #####


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####


report: "report/workflow.rst"


##### load rules #####


include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/coverage_plots.smk"
include: "rules/macs.smk"
include: "rules/peak_analysis.smk"
