Broadpeak overlap enrichment analysis, with the following fields:

  * ``qSample`` and ``tSample`` are the names of queryPeak and targetPeak files
  * ``qLen`` and ``tLen`` are the number of peaks in queryPeak and targetPeak
  * ``N_OL`` is the number of overlaps between queryPeak and targetPeak
  * ``pvalue`` the empirical pvalue of the overlap based on
    {{ snakemake.config["peak_overlap_enrichment"]["nshuffle"] }} random
    shuffle iterations
  * ``p.adjust`` the pvalue adjust for multiple comparisons using the 
    {{ snakemake.config["peak_overlap_enrichment"]["padjust_method"] }} method
