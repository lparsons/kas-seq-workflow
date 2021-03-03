rule bam_coverage:
    input:
        bam="results/dedup/{sample}-{unit}.bam",
        bai="results/dedup/{sample}-{unit}.bam.bai",
    output:
        "results/bam_coverage/{sample}-{unit}_rpgc_normalized_coverage.bw",
    log:
        "logs/bam_coverage/{sample}-{unit}.log",
    params:
        effectiveGenomeSize=config["ref"]["effective_genome_size"],
    threads: 8
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        bamCoverage \
                -b {input.bam:q} \
                -o {output:q} \
                --normalizeUsing RPGC \
                --effectiveGenomeSize {params.effectiveGenomeSize} \
                > {log:q} 2>&1
        """


rule compute_matrix:
    input:
        regions=config["ref"]["annotation"],
        score=expand(
            "results/bam_coverage/{unit.sample}-{unit.unit}_rpgc_normalized_coverage.bw",
            unit=units.itertuples(),
        ),
    output:
        matrix_gz="results/matrix_files/matrix.gz",
        #matrix_gz=report("results/matrix_files/matrix.gz",
        #        caption="../report/matrix.rst",
        #        category="Normalized KAS signal across transcripts"),
        matrix_tab="results/matrix_files/matrix.tab",
        matrix_bed="results/matrix_files/matrix.bed",
    log:
        "logs/deeptools/compute_matrix.log",
    params:
        # optional parameters
        extra="",
    threads: 16
    resources:
        mem="20G",
    conda:
        "../envs/deeptools.yml"
    shell:
        """
        computeMatrix \
            scale-regions \
            --regionsFileName {input.regions:q} \
            --upstream 3000 \
            --downstream 3000 \
            --regionBodyLength 6000 \
            --scoreFileName {input.score:q} \
            --outFileName {output.matrix_gz:q} \
            --outFileNameMatrix {output.matrix_tab:q} \
            --outFileSortedRegions {output.matrix_bed:q} \
            --numberOfProcessors {threads} \
            {params.extra} \
            > {log:q} 2>&1
        """


rule plot_heatmap:
    input:
        matrix_gz="results/matrix_files/matrix.gz",
    output:
        heatmap_img=report(
            "results/plot_heatmap/heatmap.png",
            caption="../report/heatmap.rst",
            category="Normalized KAS signal across transcripts",
        ),
        # optional output files
        regions="results/plot_heatmap/heatmap_regions.bed",
        heatmap_matrix="results/plot_heatmap/heatmap_matrix.tab",
    log:
        "logs/deeptools/heatmap.log",
    params:
        # optional parameters
        "--plotType=fill ",
    resources:
        mem="100G",
    wrapper:
        "0.72.0/bio/deeptools/plotheatmap"


rule plot_profile:
    input:
        matrix_gz="results/matrix_files/matrix.gz",
    output:
        plot_img=report(
            "results/plot_profile/profile.png",
            caption="../report/profile.rst",
            category="Normalized KAS signal across transcripts",
        ),
        # optional output files
        regions="results/plot_profile/regions.bed",
        data="results/plot_profile/data.tab",
    log:
        "logs/deeptools/plot_profile.log",
    params:
        # optional parameters
        "",
    resources:
        mem="50G",
    wrapper:
        "0.72.0/bio/deeptools/plotprofile"
