def get_fastq_group(wildcards):
    if (wildcards.group == 1):
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    else:
        return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna()

def get_fastqc(wildcards):
    for unit in units.itertuples():
        if is_single_end(unit.sample, unit.unit):
            return("qc/fastqc/{sample}-{unit}_fastqc.zip".format(sample=unit.sample, unit=unit.unit))
        else:
            return(("qc/fastqc/{sample}-{unit}.1_fastqc.zip".format(sample=unit.sample, unit=unit.unit),
                    "qc/fastqc/{sample}-{unit}.2_fastqc.zip".format(sample=unit.sample, unit=unit.unit)))

rule fastqc:
    input:
        unpack(get_fastq_group)
    output:
        html="qc/fastqc/{sample}-{unit}.{group}.html",
        zip="qc/fastqc/{sample}-{unit}.{group}_fastqc.zip"
    params:
        config["params"]["fastqc"]
    threads:
        1
    wrapper:
        "0.67.0/bio/fastqc"

rule multiqc:
    input:
        get_fastqc,
        #expand("star/{unit.sample}-{unit.unit}/Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        #expand("qc/rseqc/{unit.sample}-{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        #expand("logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples())
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.67.0/bio/multiqc"
