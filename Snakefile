subworkflow jacob_workflow:
    workdir: "jacob/diversity-data"
    snakefile: "16S.snakefile"
    configfile: "jacob/diversity-data/config.yaml"

rule all:
    input:
        expand("{x}/results.txt", x=["rossen", "moayyedi", "paramsothy", "costello", "jacob"]),
        expand("{x}/plot.pdf", x=["jacob"])

rule simple:
    wildcard_constraints: x="(rossen|moayyedi)"
    output: "{x}/results.txt"
    input: "utils.R", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

rule simple_with_data:
    wildcard_constraints: x="(paramsothy|costello)"
    output: "{x}/results.txt"
    input: "utils.R", "{x}/data.tsv", script="{x}/analyze.R"
    shell: "cd {wildcards.x} && ./analyze.R > results.txt"

rule jacob:
    output: text="jacob/results.txt", plot="jacob/plot.pdf"
    input:
        "utils.R", "jacob/patient-data.tsv",
        jacob_workflow("alpha-diversity.tsv"),
        script="jacob/analyze.R"
    shell: "cd jacob && ./analyze.R > results.txt"
