# ChIP-seq analysis pipeline

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
bowtie2-build --threads 8 -f hg38.fa hg38
index="hg38"

mkdir -p trim peaks
trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 -o  trim/ ${raw_fq}
clean_fq="trim/${sample}_trimmed.fq.gz"

# alignment
bowtie2 -p 5 -x ${index} -U ${clean_fq} | samtools sort -O bam -@ 8 -o ${sample}.bam
samtools index ${sample}.bam

# remove PCR duplicate
sambamba markdup -r -p -t 8 ${sample}.bam ${sample}.rmdup.bam
samtools index ${sample}.rmdup.bam
samtools view -h -q 30 ${sample}.rmdup.bam | grep -v chrM | samtools sort -O bam -@ 8 -o - > ${sample}.rmdup.rmchrM.bam
samtools index ${sample}.rmdup.rmchrM.bam
bedtools bamtobed -i ${sample}.rmdup.rmchrM.bam  > ${sample}.bed

# filter blacklist
bedtools intersect -v -a ${sample}.rmdup.rmchrM.bam -b hg38.blacklist.bed | samtools sort -O bam -@ 5 -o - > ${sample}.blacklist_filtered.last.bam
samtools index ${sample}.blacklist_filtered.last.bam
bedtools bamtobed -i ${sample}.blacklist_filtered.last.bam > ${sample}.last.bed

# call peaks
# Execute Processing for IP and Input/Control
# process_chip_sample <sample_prefix> <raw_fastq_path>
process_chip_sample "${chip_name}" "${chip_raw_fq}"
process_chip_sample "${control_name}" "${control_raw_fq}"
macs3 callpeak -t ${chip_name}.blacklist_filtered.last.bam -c ${control_name}.blacklist_filtered.last.bam -f BAM -n ${name} -g hs --outdir peaks/ -q 0.05

# generate bigwig file for visualization
bamCoverage --normalizeUsing RPKM -b ${sample}.blacklist_filtered.last.bam -o ${sample}.chip.bw
