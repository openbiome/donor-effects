rule all:
    input: expand("{x}/results.txt", x=["rossen", "moayyedi", "paramsothy", "costello"])

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
