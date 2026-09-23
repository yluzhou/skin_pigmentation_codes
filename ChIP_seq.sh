# ChIP-seq analysis pipeline

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
bowtie2-build --threads 8 -f hg38.fa hg38
index="hg38"

# input data
sample="sample_name"
raw_fq1="${sample}_1.fq.gz"  # Single-end input file (or R1 for paired-end)
raw_fq2="${sample}_2.fq.gz"  # R2 for paired-end

mkdir -p trim align

# Determine whether it is single-end or paired-end
if [ -n "$raw_fq2" ] && [ -f "$raw_fq2" ]; then
    echo "=== Paired-End data，start processing... ==="
    
    # Paired-end data quality control
    trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired -o trim/ ${raw_fq1} ${raw_fq2}
    
    clean_fq1="trim/${sample}_1_val_1.fq.gz"
    clean_fq2="trim/${sample}_2_val_2.fq.gz"
    
    # Paired-end reads alignment
    bowtie2 -p 5 --very-sensitive -X 2000 -x ${index} -1 ${clean_fq1} -2 ${clean_fq2} | samtools sort -O bam -@ 5 -o align/${sample}.bam

else
    echo "=== Single-End data，start processing... ==="
    
    # Single-end data quality control
    trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 -o trim/ ${raw_fq1}
    
    clean_fq="trim/${sample}_trimmed.fq.gz"
    
    # Single-end reads alignment
    bowtie2 -p 5 -x ${index} -U ${clean_fq} | samtools sort -O bam -@ 5 -o align/${sample}.bam
fi

samtools index align/${sample}.bam



# remove PCR duplicate
sambamba markdup -r -p -t 8 align/${sample}.bam align/${sample}.rmdup.bam
samtools index align/${sample}.rmdup.bam
samtools view -h -q 30 align/${sample}.rmdup.bam | grep -v chrM | samtools sort -O bam -@ 8 -o - > align/${sample}.rmdup.rmchrM.bam
samtools index align/${sample}.rmdup.rmchrM.bam
bedtools bamtobed -i align/${sample}.rmdup.rmchrM.bam  > align/${sample}.bed

# filter blacklist
bedtools intersect -v -a align/${sample}.rmdup.rmchrM.bam -b hg38.blacklist.bed | samtools sort -O bam -@ 5 -o - > align/${sample}.blacklist_filtered.last.bam
samtools index align/${sample}.blacklist_filtered.last.bam
bedtools bamtobed -i align/${sample}.blacklist_filtered.last.bam > align/${sample}.last.bed

# call peaks
# Execute Processing for IP and Input/Control
# process_chip_sample <sample_prefix> <raw_fastq_path>
process_chip_sample "${chip_name}" "${chip_raw_fq}"
process_chip_sample "${control_name}" "${control_raw_fq}"
macs3 callpeak -t ${chip_name}.blacklist_filtered.last.bam -c ${control_name}.blacklist_filtered.last.bam -f BAMPE -n ${name} -g hs --outdir peaks/ -q 0.05  # call peaks for paired-end reads
macs3 callpeak -t ${chip_name}.blacklist_filtered.last.bam -c ${control_name}.blacklist_filtered.last.bam -f BAM -n ${name} -g hs --outdir peaks/ -q 0.05  # call peaks for single-end reads

# generate bigwig file for visualization
bamCoverage --normalizeUsing RPKM -b ${sample}.blacklist_filtered.last.bam -o ${sample}.chip.bw
