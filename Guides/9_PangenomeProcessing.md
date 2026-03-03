# Supplementary Note 9: Pangenome Processing

## Alignment using GraphAligner
`GraphAligner -g hprc-v1.0-minigraph-chm13.gfa -f reads_1.fastq reads_2.fastq -a sample_ga.gam -x vg -t 20`
`GraphAligner -g hprc-v1.0-minigraph-grch38.gfa -f reads_1.fastq reads_2.fastq -a sample_ga.gam -x vg -t 20`

## Conversion to sam/bam/cram
`vg surject -x hprc-v1.0-minigraph-chm13.xg -b sample_ga.gam > sample_ga.bam`
`vg surject -x hprc-v1.0-minigraph-grch38.xg -b sample_ga.gam > sample_ga.bam`
