# RNA-seq of K562 with ActD treatment

# Requirements

    bowtie2
    STAR
    samtools

# Call SNPs

zcat results/gatk/haplotype/snps/snps.vcf.gz | grep -v '#' | awk -v OFS='\t' '{print $1,$2-1,$2 > "snps/"$1".bed"}'
for f in snps/*.bed; do bgzip $f; done

