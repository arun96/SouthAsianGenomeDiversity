# Supplementary Guide 5: RNA-Seq alignments

## STAR alignment of MAGE data to assembled contigs:

### Index the target sequences (combined set of all SAS MEGAHIT post-filtering assemblies)
`STAR --runThreadN 32 --runMode genomeGenerate --limitGenomeGenerateRAM 125000000000 --genomeSAindexNbases 13 --genomeDir RNASeq/indexes/megahit_unaligned_1kb_tagged/ --genomeFastaFiles RNASeq/combined_seqs/megahit_unaligned_1kb_tagged.fa`

### Align
`STAR --runThreadN 32 --genomeDir RNASeq/indexes/megahit_unaligned_1kb_tagged/ --readFilesIn SA_R1.fastq.gz SA_R2.fastq.gz --outFileNamePrefix RNASeq/alignments/ --readFilesCommand gunzip -c`
