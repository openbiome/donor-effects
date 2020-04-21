# Paired end read processing

include: "16S_common.snakefile"

rule paired_denoise:
    output:
        table = "table.qza",
        seqs = "rep-seqs.qza",
        stats = "denoise-stats.qza"
    input:
        "filter.qza"
    params:
        trim_length = 253,
        min_reads = 1,
        jobs_to_start = 2
    shell:
        "qiime deblur denoise-16S"
        " --i-demultiplexed-seqs {input}"
        " --p-trim-length {params.trim_length}"
        " --p-min-reads {params.min_reads}"
        " --p-jobs-to-start {params.jobs_to_start}"
        " --p-sample-stats"
        " --o-table {output.table}"
        " --o-representative-sequences {output.seqs}"
        " --o-stats {output.stats}"

rule paired_filter:
    output:
        filter = "filter.qza",
        stats = "filter-stats.qza"
    input:
        "join.qza"
    shell:
        "qiime quality-filter q-score-joined"
        " --i-demux {input}"
        " --o-filtered-sequences {output.filter}"
        " --o-filter-stats {output.stats}"

rule join:
    output:
        "join.qza"
    input:
        "demux.qza"
    shell:
        "qiime vsearch join-pairs"
        " --i-demultiplexed-seqs {input}"
        " --o-joined-sequences {output}"

rule paired_demultiplex:
    output:
        "demux.qza"
    input:
        fastq_files,
        manifest = "manifest.csv"
    shell:
        "qiime tools import"
        " --type 'SampleData[PairedEndSequencesWithQuality']"
        " --input-path {input.manifest}"
        " --input-format PairedEndFastqManifestPhred33"
        " --output-path {output}"
