#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
species_list = ["hg", "gy"]
types = ["all", "rmdup"]
indir = "results/mapping/filtered"
outdir = "results/expression"

rule all:
    input:
        expand(outdir + "/fpkm/{sample}.{species}.{t}.tsv", sample=samples, species=species_list, t=types),
        expand(outdir + "/feature_count/{sample}.{species}.{t}.tsv", sample=samples, species=species_list, t=types),

def get_bam(wildcards):
    if wildcards.t == "all":
        return "results/mapping/filtered/%s.%s.bam" % (wildcards.sample, wildcards.species)
    else:
        return "results/mapping/rmdup/%s.%s.bam" % (wildcards.sample, wildcards.species)

def get_bed(species):
    if species == "hg":
        return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.transcripts.bed.gz"
    else:
        return "/date/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.transcripts.bed.gz"

def get_tsv(species):
    if species == "hg":
        return "/date/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.tsv"
    else:
        return "/date/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.annotation.tsv"

rule calculate_fpkm:
    input:
        bam = lambda wildcards: get_bam(wildcards),
        bed = lambda wildcards: get_bed(wildcards.species),
        tsv = lambda wildcards: get_tsv(wildcards.species)
    output:
        tmp = temp(outdir + "/fpkm/{sample}.{species}.{t}.TMP.tsv"),
        tsv = outdir + "/fpkm/{sample}.{species}.{t}.tsv"
    log:
        outdir + "/fpkm/{sample}.{species}.{t}.log"
    threads:
        8
    shell:
        """(
        ./scripts/expression/calculate_fpkm.stranded.v2.py {input.bam} {input.bed} PE {threads} {output.tmp}
        ./scripts/expression/annotate_fpkm.py {output.tmp} {input.tsv} {output.tsv} ) &> {log}
        """

def get_gtf(species):
    if species == "hg":
        return "/home/chenzonggui/species/homo_sapiens/GRCh38.p13/gencode.v39.annotation.sorted.gtf"
    else:
        return "/home/chenzonggui/species/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.108.gtf"

rule feature_count:
    input:
        bam = lambda wildcards: get_bam(wildcards),
        gtf = lambda wildcards: get_gtf(wildcards.species)
    output:
        tsv = outdir + "/feature_count/{sample}.{species}.{t}.tsv"
    log:
        outdir + "/feature_count/{sample}.{species}.{t}.log"
    threads:
        8
    shell:
        """
        featureCounts -T {threads} -s 2 -p -B -a {input.gtf} -o {output.tsv} {input.bam} &> {log}
        """