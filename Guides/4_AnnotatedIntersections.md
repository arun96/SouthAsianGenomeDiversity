# Supplementary Guide 4: Intersection with Annotated Elements

## Bedtools intersection with annotated elements:

### Bam to Bed conversion
`samtools view -Sb megahit_unaligned_1kb_tagged_mapped.sam | bedtools bamtobed -i - | bedtools sort -i - > megahit_unaligned_1kb_tagged_mapped.sorted.bed`

### Intersect with annotations file - keep details of intersection
`bedtools intersect -a megahit_unaligned_1kb_tagged_mapped.sorted.bed -b chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz -wb > megahit_unaligned_1kb_tagged_mapped.intersect.bed`

### No details, simple .bed file
`bedtools intersect -a megahit_unaligned_1kb_tagged_mapped.sorted.bed -b chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz > megahit_unaligned_1kb_tagged_mapped.intersect_simple.bed`

## Bedtools intersection with GWAS sites (for direct overlap):

### Intersect with annotations file - keep details of intersection
`bedtools intersect -a megahit_unaligned_1kb_tagged_mapped.sorted.bed -b chm13v2.0_GWASv1.0rsids_e100_r2022-03-08.vcf.gz -wb > megahit_unaligned_1kb_tagged_mapped.intersect_gwas.bed`

### No details, simple .bed file
`bedtools intersect -a megahit_unaligned_1kb_tagged_mapped.sorted.bed -b chm13v2.0_GWASv1.0rsids_e100_r2022-03-08.vcf.gz > megahit_unaligned_1kb_tagged_mapped.intersect_gwas_simple.bed`
