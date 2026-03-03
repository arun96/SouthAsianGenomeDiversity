# Supplementary Guide 8: All vs all alignment to identify shared sequences

## Self alignment with no limit on the number of secondary alignments
### Option 1
`minimap2 -c -x asm5 -N1000 --cs -t 24 all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa` `all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa > all_sas_1kb_tagged_no_outliers_nolimit.paf`

### Option 2
`minimap2 -DP -cx asm5 -t 24 all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa > all_sas_1kb_tagged_no_outliers_ava.paf`