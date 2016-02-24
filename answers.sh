 #! /usr/bin/env bash
datasets='/vol3/home/johnsonr/data-sets'

#Question 1(answer: 1079): Use BEDtools intersect to identify the size of the largest overlap between CTCF and H3K4me3 locations.
chr22bed="$datasets/bed/encode.tfbs.chr22.bed.gz"
h3k4bed="$datasets/bed/encode.h3k4me3.hela.chr22.bed"
answer_1=$(zcat $chr22bed \
| sort -k4 \
| awk '($4 == "CTCF")' \
| bedtools intersect -a stdin -b $h3k4bed -wo \
| awk '{print $NF}' \
| sort -nr \
| head -n1 )

echo "answer-1: $answer_1"

#Question 2 (answer: 0.384000): Use BEDtools to calculate the GC content of nucleotides 19,000,000 to 19,000,500 on chr22 of hg19 genome build. Report the GC content as a fraction (e.g., 0.50).
fasta19="$datasets/fasta/hg19.chr22.fa"
answer_2=$( echo -e "chr22\t19000000\t19000500" > $datasets/fasta/chr22stsp.bd \
| bedtools nuc -fi $fasta19 -bed $datasets/fasta/chr22stsp.bd \
| cut -f5 \
| tail -n1)

echo "answer-2: $answer_2"

#Q3(answer: 850): Use BEDtools to identify the length of the CTCF ChIP-seq peak (ie interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz
# calculate signal in intervals, bedtools map, summary output stat=mean
h3k4bg="$datasets/bedtools/ctcf.hela.chr22.bg.gz"
answer_3=$(bedtools map -a $chr22bed  -b $h3k4bg -c 4 -o mean \
| sort -k4 \
| awk '($4 == "CTCF")' \
| sort -k5nr \
| awk '(NR==1) {interval = $3 - $2} END {print interval}')

echo "answer-3: $answer_3"

#Q4(PRAME): use bedtools to identify the gene promoter (defined as 1000 bp upstream of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz. Report the gene name (e.g., 'ABC123') tss file, find promoters (1000 bases upstream of tss), bedtools flank, use -s flag to pay attention to strand, compare intervals to signal using map, summary stat=median, print gene name
tss="$datasets/bed/tss.hg19.chr22.bed.gz"
genome="$datasets/genome/hg19.genome"
answer_4=$(sort $tss -k2n \
| bedtools flank -i $tss -g $genome -l 1000 -r 0 -s \
| sort -k2n \
| bedtools map -a stdin -b $h3k4bg -c 4 -o median \
| sort -k7nr \
| awk '(NR==1) {print $4}')

echo "answer-4: $answer_4"

#Q5(answer: "chr22:0-16150259"): bedtools complement (intersect -v should also work--but more complex), calculate and sort by new interval length, but report interval
hg19genes="$datasets/bed/genes.hg19.bed.gz"
answer_5=$(zcat $hg19genes \
| grep "chr22" \
| sort -k2n \
| bedtools complement -i stdin -g $genome \
| awk '($1=="chr22") {print $1, $2, $3, $3-$2}' \
| sort -k4nr \
| awk 'BEGIN {OFS=""} (NR==1) {print $1,":",$2,"-",$3}')

echo "answer-5: $answer_5"

answer_6=$(bedtools jaccard -a $chr22bed -b $h3k4bed \
| sort -k3nr \
| awk '(NR==1) {print $3}')

echo "answer-6: $answer_6"
