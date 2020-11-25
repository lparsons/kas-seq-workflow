MACS2 peaks - A tabular file which contains information about called peaks. You can open it in excel and sort/filter using excel functions. Information include:

* chromosome name
* start position of peak
* end position of peak
* length of peak region
* absolute peak summit position
* pileup height at peak summit
* ``-log10(pvalue)`` for the peak summit (*e.g.* ``pvalue =1e-10``, then this value should be ``10``)
* fold enrichment for this peak summit against random Poisson distribution with local lambda,
* ``-log10(qvalue)`` at peak summit

Coordinates in XLS are 1-based which is different from BED format. Note that since ``--broad`` is enabled, the pileup, p-value, q-value, and fold change in the XLS file will be the mean value across the entire peak region, since peak summit won't be called in broad peak calling mode.
