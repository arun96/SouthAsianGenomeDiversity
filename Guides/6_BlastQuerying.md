# Supplementary Guide 6: BLAST querying

## BLAST querying against nt db

### Unplaced contigs
`blastn -db nt -query megahit_unaligned_1kb_tagged_unplaced.fa -out megahit_unaligned_1kb_tagged_unplaced.out -num_threads 32`
### Placed contigs
`blastn -db nt -query megahit_unaligned_1kb_tagged_placed.fa -out megahit_unaligned_1kb_tagged_placed.out -num_threads 32`

### Extract the top 50 hits for each sequence
`grep 'Query=' megahit_unaligned_1kb_tagged_unplaced.out -A 50 > megahit_unaligned_1kb_tagged_unplaced_summary50.out`
`grep 'Query=' megahit_unaligned_1kb_tagged_placed.out -A 50 > megahit_unaligned_1kb_tagged_placed_summary50.out`

Post process to analyze these hits.
