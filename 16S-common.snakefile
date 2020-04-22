import csv
import hashlib
import os.path, wget
from snakemake.exceptions import IncompleteCheckpointException

# Functions -----------------------------------------------------------

def download(url, filepath, md5):
    if not os.path.isfile(filepath):
        wget.download(url, filepath)

    with open(filepath, "rb") as f:
        observed_md5 = hashlib.md5(f.read()).hexdigest()
        if not observed_md5 == md5:
            raise RuntimeError(f"Downloaded file {filepath} from url {url}"
                "has MD5 {observed_md5}, which does not match file manifest's"
                "expected MD5 {md5}")

def minimum_depth(fn):
    '''Get minimum number of counts in OTU table across samples'''
    with open(fn) as f:
        next(f) # skip first line
        reader = csv.DictReader(f, dialect='excel-tab')
        rows = [x for x in reader]

    samples = [x for x in rows[0].keys() if x != "#OTU ID"]
    sample_counts = [sum([float(row[sample]) for row in rows]) for sample in samples]
    return min(sample_counts)


# Rules ---------------------------------------------------------------

rule export_beta:
    output: "distance-matrix.tsv"
    input: "distance-matrix.qza"
    shell: "qiime tools export --input-path {input} --output-path ."

rule beta:
    output: "distance-matrix.qza"
    input: "rarefied-table.qza"
    shell:
        "qiime diversity beta"
        " --i-table {input}"
        " --o-distance-matrix {output}"
        " --p-metric braycurtis"

rule export_diversity:
    output:
        "alpha-diversity.tsv"
    input:
        "alpha-diversity.qza"
    shell:
        "qiime tools export --input-path {input} --output-path ."

rule diversity:
    output:
        "alpha-diversity.qza"
    input:
        "rarefied-table.qza"
    params:
        metric = "shannon"
    shell:
        "qiime diversity alpha"
        " --i-table {input}"
        " --p-metric {params.metric}"
        " --o-alpha-diversity {output}"

rule rarefy:
    output:
        "rarefied-table.qza"
    input:
        table = "table.qza",
        depth_file = "rarefaction-depth.txt"
    shell:
        "qiime feature-table rarefy"
        " --i-table {input.table}"
        " --p-sampling-depth $(perl -ne '/(\d+)\./; print $1;' {input.depth_file})"
        " --o-rarefied-table {output}"

rule compute_depth:
    output:
        "rarefaction-depth.txt"
    input:
        "table.tsv"
    run:
        with open(output[0], 'w') as f:
            f.write(str(minimum_depth(input[0])))

rule convert_table:
    output:
        "table.tsv"
    input:
        "feature-table.biom"
    shell:
        "biom convert --to-tsv -i {input} -o {output}"

rule export_table:
    output:
        "feature-table.biom"
    input:
        "table.qza"
    shell:
        "qiime tools export --input-path {input} --output-path ."

def fastq_files(wildcards):
    if os.path.exists("samples.csv"):
        with open("samples.csv") as f:
            reader = csv.DictReader(f)
            filepaths = [row["filepath"] for row in reader]

        return filepaths
    else:
        raise IncompleteCheckpointException("manifest", "samples.csv")

rule download:
    output:
        "{x}.fastq.gz"
    input:
        "samples.csv"
    run:
        with open(input[0]) as f:
            reader = csv.DictReader(f)
            row = [row for row in reader if row["filepath"] == output[0]]

        assert len(row) == 1
        row = row[0]

        download(row["url"], row["filepath"], row["md5"])

checkpoint manifest:
    output:
        manifest = "manifest.csv",
        samples = "samples.csv"
    input:
        filereport = "filereport.tsv",
        metadata = config["metadata"],
        script = "make-manifest.R"
    shell:
        "./{input.script} --metadata {input.metadata} --filereport {input.filereport}"
        " --manifest {output.manifest} --samples {output.samples}"

rule filereport:
    params:
        accession = config["study-accession"]
    output:
        "filereport.tsv"
    shell:
        "wget https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={params.accession}\\&result=read_run -O {output}"
