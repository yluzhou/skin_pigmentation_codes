# ATAC-seq analysis pipeline
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip -c hg38.fa.gz > hg38.fa
bowtie2-build --threads 8 -f hg38.fa hg38
index="hg38"

mkdir -p trim peaks
trim_galore --illumina -q 25 --phred33 --length 35 -e 0.1 --stringency 4 --paired -o trim/ ${fq1} ${fq2}
clean_fq1="trim/${sample}_1_val_1.fq.gz"
clean_fq2="trim/${sample}_2_val_2.fq.gz"

# alignment
bowtie2 -p 5 --very-sensitive -X 2000 -x ${index} -1 ${clean_fq1} -2 ${clean_fq2} | samtools sort -O bam -@ 8 -o - > ${sample}.bam
samtools index ${sample}.bam

# remove PCR duplicates 
sambamba markdup --overflow-list-size 600000  --tmpdir='./'  -r ${sample}.bam ${sample}.rmdup.bam
samtools index ${sample}.rmdup.bam
samtools idxstats ${sample}.rmdup.bam | cut -f 1 | grep -v -E 'chrM|chrMT|MT' | xargs samtools view -b -f 2 -q 30 -@ 5 ${sample}.rmdup.bam > ${sample}.rmdup.rmchrM.bam
samtools index ${sample}.rmdup.rmchrM.bam
bedtools bamtobed -i ${sample}.rmdup.rmchrM.bam  > ${sample}.bed

# filter blacklist
bedtools intersect -v -a ${sample}.rmdup.rmchrM.bam -b hg38.blacklist.bed | samtools sort -O bam -@ 8 -o - > ${sample}.blacklist_filtered.last.bam
samtools index ${sample}.blacklist_filtered.last.bam
bedtools bamtobed -i ${sample}.blacklist_filtered.last.bam > ${sample}.last.bed

# call peaks
macs3 callpeak -t ${sample}.blacklist_filtered.last.bam -g hs -f BAMPE -n ${sample_name} --outdir ${out}/ -q 0.05

# generate bigwig file for visualization
bamCoverage --normalizeUsing RPKM -b ${sample}.blacklist_filtered.last.bam -o ${sample}.atac.bw
