# Single-end reads

include: "16S-common.snakefile"

rule single_denoise:
    output:
        table = "table.qza",
        seqs = "rep-seqs.qza",
        stats = "denoise-stats.qza"
    input:
        "filter.qza"
    params:
        trim_length = 150,
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

rule single_filter:
    output:
        filter = "filter.qza",
        stats = "filter-stats.qza"
    input:
        "demux.qza"
    shell:
        "qiime quality-filter q-score"
        " --i-demux {input}"
        " --o-filtered-sequences {output.filter}"
        " --o-filter-stats {output.stats}"

rule single_demultiplex:
    output:
        "demux.qza"
    input:
        fastq_files,
        manifest = "manifest.csv"
    shell:
        "qiime tools import"
        " --type 'SampleData[SequencesWithQuality]'"
        " --input-path {input.manifest}"
        " --input-format SingleEndFastqManifestPhred33"
        " --output-path {output}"
