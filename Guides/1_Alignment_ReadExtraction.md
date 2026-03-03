# Supplementary Guide 1: Alignment and Read Extraction

## Alignment of input reads to chosen reference:

### Index reference
`bowtie2-build ref.fa indexes/ref --threads ${THRDS}`

### Align reads
`bowtie2 -p ${THRDS} -x indexes/ref -1 ${RDS}_1.fastq -2 ${RDS}_2.fastq > ${EXP}.sam`

### Post-processing
`samtools view -@ ${THRDS} -S -b ${EXP}.sam > ${EXP}.bam`
`samtools sort -@ ${THRDS} ${EXP}.bam -o ${EXP}.sorted.bam`


## Extraction of unaligned read pairs, and unaligned reads with a mapped mate:

### Get unaligned reads
`samtools fastq -f 12 $ALGN -1 unmapped_R1.fq -2 unmapped_R2.fq -@ ${THRDS}`

### Get first reads that are unmapped but have a mapped mate
`samtools fastq -f 68 -F 8 $ALGN > mateMapped_R1.fq -@ ${THRDS}`

### Get second reads that are unmapped but have a mapped mate
`samtools fastq -f 132 -F 8 $ALGN > mateMapped_R2.fq -@ ${THRDS}`

### Get all unaligned reads into one file
`echo "Extracting All Unaligned Reads"`
`samtools view -b -f 4 $ALGN > unmapped.bam -@ ${THRDS}`
`samtools sort unmapped.bam -o unmapped_sorted.bam -@ ${THRDS}`

### Extraction of reads based on quality value:
`samtools view -h $ALGN -q 20 -o greater_than_q20.sam -U q20.sam -@ ${THRDS}`
`samtools sort q20.sam -o q20_sorted.sam -@ ${THRDS}`
`samtools fastq q20_sorted.sam -1 q20_R1.fq -2 q20_R2.fq -@ ${THRDS}`