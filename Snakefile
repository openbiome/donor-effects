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
        expand("{x}/results.txt", x=["rossen", "moayyedi", "paramsothy", "costello", "jacob", "goyal", "kump", "nishida", "pools"]),
        expand("{x}/plot.pdf", x=["jacob", "kump"])

rule simple:
    wildcard_constraints: x="(rossen|moayyedi|nishida|pools)"
    output: "{x}/results.txt"
    input: "utils.R", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

rule simple_with_data:
    wildcard_constraints: x="(paramsothy|costello)"
    output: "{x}/results.txt"
    input: "utils.R", "{x}/data.tsv", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

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
