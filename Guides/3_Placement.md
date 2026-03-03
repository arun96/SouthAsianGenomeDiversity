# Supplementary Guide 3: Placement

## Minimap2 alignment of contigs to the chosen reference:

### Alignment in SAM format
`minimap2 -a -t 16 ref.fa megahit_unaligned_1kb_tagged.fa > megahit_unaligned_1kb_tagged_mapped.sam`

### Alignment in PAF format
`minimap2 -c -t 16 ref.fa megahit_unaligned_1kb_tagged.fa > megahit_unaligned_1kb_tagged_mapped.paf`

In addition to this, we use the mate-pair linking approach taken in the African Pangenome effort (Sherman et al, 2019).
