#!/usr/bin/env runsnakemake
include: "0_SnakeCommon.smk"
indir = "results/mapping/rmdup"
outdir = "results/introns"

rule all:
    input:
        expand(outdir + "/counts/{sample}.tsv", sample=samples),

rule quantify_introns:
    input:
        bed = outdir + "/introns.bed.gz",
        bam = indir + "/{sample}.human.bam"
    output:
        txt = outdir + "/counts/{sample}.tsv"
    log:
        outdir + "/counts/{sample}.log"
    threads:
        8
    shell:
        """
        ./scripts/quantify_introns.py {input.bed} {input.bam} {threads} {output.txt} &> {log}
        """