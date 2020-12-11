rule peak_annotation:
    input:
        peak_files=expand("results/macs2/{unit.sample}-{unit.unit}_peaks.broadPeak", unit=units.itertuples()),
        annotation_file=config["ref"]["annotation"],
    output:
        annotation_distribution_plot=report(
                "results/peak_analysis/broadpeakAnnotationDistributionPlot.pdf",
                caption="../report/broadpeakAnnotationDistributionPlot.rst",
                category="KAS-Seq Peak Analysis"),
        peak_annotation_list_rdata=report(
                "results/peak_analysis/broadPeakAnnoList.Rdata",
                caption="../report/broadpeakAnnoList.rst",
                category="KAS-Seq Peak Analysis"),
        annotated_peak_files=report(
                expand("results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.tsv.gz", unit=units.itertuples()),
                caption="../report/broadpeak_annotation.rst",
                category="KAS-Seq Peak Analysis"),
        annotated_peak_summary=report(
                expand("results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.summary.txt", unit=units.itertuples()),
                caption="../report/broadpeak_annotation_summary.rst",
                category="KAS-Seq Peak Analysis"),
    params:
        names=expand("{unit.sample}-{unit.unit}", unit=units.itertuples()),
        output_dir="results/peak_analysis",
    resources:
        mem="16G"
    log:
        "logs/peak_analysis/broadpeak_annotation.log"
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/peak_annotation.R"


rule peak_profile:
    input:
        peak_files=expand("results/macs2/{unit.sample}-{unit.unit}_peaks.broadPeak", unit=units.itertuples()),
        annotation_file=config["ref"]["annotation"],
    output:
        tag_profile_plot=report(
                "results/peak_analysis/broadPeakAvgProfilePlot.pdf",
                caption="../report/broadpeakAvgProfilePlot.rst",
                category="KAS-Seq Peak Analysis"),
        tag_heatmap_plot=report(
                "results/peak_analysis/broadPeakHeatmapPlot.pdf",
                caption="../report/broadpeakHeatmapPlot.rst",
                category="KAS-Seq Peak Analysis"),
        tag_matrix_list="results/peak_analysis/broadPeakMatrixList.Rdata",
    params:
        names=expand("{unit.sample}-{unit.unit}", unit=units.itertuples()),
    resources:
        mem="16G"
    log:
        "logs/peak_analysis/broadpeak_profile.log"
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/peak_profile.R"

