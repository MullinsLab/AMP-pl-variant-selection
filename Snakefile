configfile: "config.yaml"
BLAST = "/usr/local/bin" #configure path to blast executable
DATASETS = [d for d in config for s in config[d]]
SAMPLES = [s for d in config for s in config[d]]

##Possible to run porpid-postproc as a subworkflow
##Simply add porpid_postproc() calls to linked input
#subworkflow porpid_postproc:
#    workdir:
#        "../porpid-postproc"
#    snakefile:
#        "../porpid-postproc/Snakefile"
#    configfile:
#        "../porpid-postproc/config.yaml"

def variant_params(wc):
    if "REN" in wc["sample"]:
        ref = "panels/HXB2_5970-9012.fasta"
        start = 256
        end = 2826
    elif "GP" in wc["sample"]:
        ref = "panels/HXB2_757-3300.fasta"
        start = 33
        end = 2322
    else: #defauts to REN
        ref = "panels/HXB2_5970-9012.fasta"
        start = 256
        end = 2826
    blast = BLAST
    return {"ref": ref, "start": start, "end": end, "blast": blast}

rule all:
    input:
        expand("variant_analysis/{dataset}/{sample}/{sample}_log.txt", zip, dataset = DATASETS, sample = SAMPLES)

rule analyze_variants:
    input:
        "postproc/{dataset}/{sample}/{sample}.fasta" #add porpid_postproc() here if running subworkflow
    params:
        p = variant_params
    output:
        "variant_analysis/{dataset}/{sample}/{sample}_log.txt"
    shell:
        "perl scripts/runProcess.pl -sid {wildcards.sample} -in {input} -envref {params[p][ref]} -rs {params[p][start]} -re {params[p][end]} -outpath variant_analysis/{wildcards.dataset}/{wildcards.sample}/ -blast {params[p][blast]}"
