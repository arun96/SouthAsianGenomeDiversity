# Supplementary Guide 7: Long Read Assembly and Validation

## Long read assembly using Flye:
### Convert unaligned read .bam to .fasta
`samtools -@ 42 fasta sample.ONT.unaligned.bam > sample_ONT.fasta`

### Flye assembly
`mkdir flye_ONT_raw`
`flye --nano-raw sample_ONT.fasta --out-dir ./flye_ONT_raw --threads 42 -g 3.1g --debug | tee flye.log`

## Alignment of unaligned reads to long read assembly:
### Index
`bowtie2-build sample_ONT.fasta indexes/sample_ONT --threads 36`

### Align
`(bowtie2 -p 24 -x indexes/sample_ONT -1 unmapped_R1.fq -2 unmapped_R2.fq -S unmapped.sam) 2> unmapped_ONT.log`

## Alignment of assembled contigs to long read assembly:
#### Alignment in sam/paf format
`minimap2 -a -t 32 ONT_assembly.fasta sample_contigs_1kb.fa > unaligned_1kb_contigs.sam `
`minimap2 -t 32 ONT_assembly.fasta sample_contigs_1kb.fa > unaligned_1kb_contigs.paf`

#### Counting total contigs
`grep '>' sample_contigs_1kb.fa | wc -l`

#### Counting alignments with alignment length > 1Kb
`awk '$11>1000' unaligned_1kb_contigs.paf | cut -f 1 | sort | uniq | wc -l`
