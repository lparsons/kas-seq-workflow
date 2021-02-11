def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_fastq(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if is_single_end(wildcards.sample, wildcards.unit):
        return {"fq1": f"{u.fq1}"}
    else:
        return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        if is_single_end(wildcards.sample, wildcards.unit):
            return [f"{u.fq1}"]
        else:
            return [f"{u.fq1}", f"{u.fq2}"]

    else:
        # yes trimming, use trimmed data
        if not is_single_end(wildcards.sample, wildcards.unit):
            # paired-end sample
            return expand(
                "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
                group=[1, 2],
                sample=wildcards.sample,
                unit=wildcards.unit,
            )
        # single end sample
        return [
            "results/trimmed/{sample}-{unit}.fastq.gz".format(
                sample=wildcards.sample, unit=wildcards.unit
            )
        ]


def get_fastq_group(wildcards):
    fastq_files = get_fq(wildcards)
    return fastq_files[int(wildcards.group) - 1]


def get_fastqc(wildcards):
    fastqc_files = list()
    for unit in units.itertuples():
        if is_single_end(unit.sample, unit.unit):
            fastqc_files.append(
                "qc/fastqc/{sample}-{unit}_fastqc.zip".format(
                    sample=unit.sample, unit=unit.unit
                )
            )
        else:
            fastqc_files.extend(
                [
                    "qc/fastqc/{sample}-{unit}.1_fastqc.zip".format(
                        sample=unit.sample, unit=unit.unit
                    ),
                    "qc/fastqc/{sample}-{unit}.2_fastqc.zip".format(
                        sample=unit.sample, unit=unit.unit
                    ),
                ]
            )
    return fastqc_files


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
