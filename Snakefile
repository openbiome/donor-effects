STUDIES=["rossen", "moayyedi", "paramsothy", "costello", "jacob", "goyal", "kump", "nishida", "pools"]

subworkflow jacob_workflow:
    workdir: "jacob/diversity-data"
    snakefile: "16S.snakefile"
    configfile: "jacob/diversity-data/config.yaml"

subworkflow goyal_workflow:
    workdir: "goyal/diversity-data"
    snakefile: "16S.snakefile"
    configfile: "goyal/diversity-data/config.yaml"

subworkflow kump_workflow:
    workdir: "kump/diversity-data"
    snakefile: "16S.snakefile"
    configfile: "kump/diversity-data/config.yaml"

rule all:
    input:
        expand("{x}/results.txt", x=STUDIES),
        expand("{x}/plot.pdf", x=["jacob", "kump"])

# Three studies have their data embedded in the scripts
rule simple:
    wildcard_constraints: x="(rossen|moayyedi|nishida)"
    output: "{x}/results.txt"
    input: "utils.R", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

# Another 2 have separate data files
rule simple_with_data:
    wildcard_constraints: x="(paramsothy|costello)"
    output: "{x}/results.txt"
    input: "utils.R", "{x}/patient-data.tsv", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

# Pool meta-analysis depends on other studies
rule pools:
    output: "pools/results.txt"
    input:
        "utils.R", "pools/analyze.R",
        expand("{x}/patient-data.tsv", x=["paramsothy", "costello", "jacob"])
    shell: "cd pools && ./analyze.R > results.txt"

# The three microbiome
rule jacob:
    output: "jacob/results.txt", "jacob/plot.pdf"
    input:
        "utils.R", "jacob/patient-data.tsv",
        jacob_workflow("alpha-diversity.tsv"),
        script="jacob/analyze.R"
    shell: "cd jacob && ./analyze.R > results.txt"

rule goyal:
    output: "goyal/results.txt"
    input:
        "utils.R", "goyal/patient-data.tsv",
        goyal_workflow("alpha-diversity.tsv"),
        script="goyal/analyze.R"
    shell: "cd goyal && ./analyze.R > results.txt"

rule kump:
    output: "kump/results.txt", "kump/plot.pdf"
    input:
        "utils.R", "kump/kump2018.metadata.tsv",
        kump_workflow("alpha-diversity.tsv"),
        script="kump/analyze.R"
    shell: "cd kump && ./analyze.R > results.txt"
