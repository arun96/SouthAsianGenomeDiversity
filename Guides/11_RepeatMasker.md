# Supplementary Guide 11: Running RepeatMasker

## Default run
`RepeatMasker --species Human all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa`

## Run without masking low complexity DNA or repeats
`RepeatMasker --species Human -nolow all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa`

Key output files are in all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa.masked, all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa.out and all_sas_megahit_unaligned_1kb_tagged_no_outliers.fa.tbl
