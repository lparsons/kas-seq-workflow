rule samtools_faindex:
    input:
        "{sample}.fa",
    output:
        "{sample}.fa.fai",
    params:
        "",  # optional params string
    wrapper:
        "0.68.0/bio/samtools/faidx"


rule make_txdb:
    input:
        txdb_input=config["ref"]["annotation"],
        chrom_info=config["ref"]["fasta"] + ".fai",
    output:
        "results/peak_analysis/txdb.db",
    log:
        "logs/peak_analysis/make_txdb.log",
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/save_txdb.R"


rule peak_annotation:
    input:
        peak_files=expand(
            "results/macs2/{unit.sample}-{unit.unit}_peaks.broadPeak",
            unit=units.itertuples(),
        ),
        txdb_file="results/peak_analysis/txdb.db",
    output:
        annotation_distribution_plot=report(
            "results/peak_analysis/broadpeakAnnotationDistributionPlot.pdf",
            caption="../report/broadpeakAnnotationDistributionPlot.rst",
            category="KAS-Seq Peak Analysis",
        ),
        peak_annotation_list_rdata=report(
            "results/peak_analysis/broadPeakAnnoList.Rdata",
            caption="../report/broadpeakAnnoList.rst",
            category="KAS-Seq Peak Analysis",
        ),
        annotated_peak_files=report(
            expand(
                "results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.tsv.gz",
                unit=units.itertuples(),
            ),
            caption="../report/broadpeak_annotation.rst",
            category="KAS-Seq Peak Analysis",
        ),
        annotated_peak_summary=report(
            expand(
                "results/peak_analysis/{unit.sample}-{unit.unit}_peaks.broadPeak.annotated.summary.txt",
                unit=units.itertuples(),
            ),
            caption="../report/broadpeak_annotation_summary.rst",
            category="KAS-Seq Peak Analysis",
        ),
        venn_diagram=report(
            "results/peak_analysis/broadPeakAnnotationVennDiagram.pdf",
            caption="../report/broadpeak_annotation_venn_diagram.rst",
            category="KAS-Seq Peak Analysis",
        ),
    params:
        names=expand("{unit.sample}-{unit.unit}", unit=units.itertuples()),
        output_dir="results/peak_analysis",
    resources:
        mem="16G",
    log:
        "logs/peak_analysis/broadpeak_annotation.log",
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/peak_annotation.R"


rule peak_profile:
    input:
        peak_files=expand(
            "results/macs2/{unit.sample}-{unit.unit}_peaks.broadPeak",
            unit=units.itertuples(),
        ),
        txdb_file="results/peak_analysis/txdb.db",
    output:
        tag_profile_plot=report(
            "results/peak_analysis/broadPeakAvgProfilePlot.pdf",
            caption="../report/broadpeakAvgProfilePlot.rst",
            category="KAS-Seq Peak Analysis",
        ),
        tag_heatmap_plot=report(
            "results/peak_analysis/broadPeakHeatmapPlot.pdf",
            caption="../report/broadpeakHeatmapPlot.rst",
            category="KAS-Seq Peak Analysis",
        ),
        tag_matrix_list="results/peak_analysis/broadPeakMatrixList.Rdata",
    params:
        names=expand("{unit.sample}-{unit.unit}", unit=units.itertuples()),
    resources:
        mem="16G",
    log:
        "logs/peak_analysis/broadpeak_profile.log",
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/peak_profile.R"


def get_target_broadpeaks(wildcards):
    query_select = (units["sample"] == wildcards.sample) & (
        units["unit"] == wildcards.unit
    )
    query_unit = units[query_select]
    target_units = units[~query_select]
    target_peak_files = expand(
        "results/macs2/{unit.sample}-{unit.unit}_peaks.broadPeak",
        unit=target_units.itertuples(),
    )
    return target_peak_files


rule peak_overlap_enrichment:
    input:
        query_peak_file="results/macs2/{sample}-{unit}_peaks.broadPeak",
        target_peak_files=get_target_broadpeaks,
        txdb_file="results/peak_analysis/txdb.db",
    output:
        output_tsv_file=report(
            "results/peak_analysis/{sample}-{unit}_broadpeak_overlap_enrichment.tsv",
            caption="../report/broadpeak_overlap_enrichment.rst",
            category="KAS-Seq Peak Analysis",
        ),
    params:
        nshuffle=config["peak_overlap_enrichment"]["nshuffle"],
        p_adjust_method=config["peak_overlap_enrichment"]["p_adjust_method"],
    threads: 1
    resources:
        mem="16G",
    log:
        "logs/peak_analysis/broadpeak_overlap_enrichment_{sample}-{unit}.log",
    conda:
        "../envs/chipseeker.yml"
    script:
        "../scripts/peak_overlap_enrichment.R"
