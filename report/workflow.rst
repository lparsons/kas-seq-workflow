This workflow performs a KAS-Seq analysis on single- or paired-end KAS-Seq data (`Wu *et al.* 2020. Kethoxal-assisted single-stranded DNA sequencing captures global transcription dynamics and enhancer activity in situ <https://doi.org/10.1038/s41592-020-0797-9>`_).

The main steps of the workflow are:

1.  Adapter removal with `Cutadapt <http://cutadapt.readthedocs.io/>`_

2.  Mapping using `BWA-MEM <http://bio-bwa.sourceforge.net/bwa.shtml>`_

3.  Deduplication with `Picard MarkDuplicates <https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates>`_

4.  `Deeptools <https://deeptools.readthedocs.io>`_ `computMatrix <https://deeptools.readthedocs.io/en/latest/content/tools/computeMatrix.html>`_,
    `plotHeatMap <https://deeptools.readthedocs.io/en/latest/content/tools/plotHeatmap.html>`_,
    and `plotFingerprint <https://deeptools.readthedocs.io/en/latest/content/tools/plotProfile.html>`_

5.  Peak calling using `MACS2 <https://github.com/macs3-project/MACS>`_

