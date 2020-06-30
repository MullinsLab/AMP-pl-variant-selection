from glob import glob
from os.path import basename
import re

BLAST = "/Users/appankow/opt/anaconda3/envs/porpid-postproc/bin/" #configure path to blast executable
FILES = glob("data/*.fasta")
SAMPLES = [re.match("(.+).fasta",basename(f))[1] for f in FILES]

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
        expand("variant_analysis/{sample}",sample = SAMPLES)

rule analyze_variants:
    input:
        "data/{sample}.fasta" #add porpid_postproc() here if running subworkflow
    params:
        p = variant_params
    output:
        directory("variant_analysis/{sample}")
    shell:
        "perl scripts/runProcess.pl -sid {wildcards.sample} -in {input} -envref {params[p][ref]} -rs {params[p][start]} -re {params[p][end]} -outpath variant_analysis/{wildcards.sample}/ -blast {params[p][blast]}"

#rule select_variants:
#    input:
#        "variant_analysis/{sample}",
#        "data/{sample}/{sample}.fasta.mafft.fa"
#    params:
#        SID = lambda wc: wc.sample
#    output:
#        "select_variants/{sample}_select_var.fasta"
#    script:
#        "scripts/select-variants.jl"
