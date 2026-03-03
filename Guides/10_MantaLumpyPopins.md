# Supplementary Guide 10: Running Manta, Lumpy and PopIns2

## Manta

### prepare the analysis job
```sh
configManta.py \
  --bam sample~{bam_ext} \
  --referenceFasta ~{reference_fasta} \
  --runDir . &&

./runWorkflow.py \
  --mode local \
  --jobs ~{num_jobs} \
  --memGb $((~{num_jobs} * 2))
```

### inversion conversion, then compression and index
```sh
python2 /usr/local/bin/manta/libexec/convertInversion.py \
  /usr/local/bin/samtools \
  ~{reference_fasta} \
  results/variants/diploidSV.vcf.gz \
  | bcftools reheader -s <(echo "~{sample_id}") \
  > diploidSV.vcf

bgzip -c diploidSV.vcf > ~{sample_id}.manta.vcf.gz
tabix -p vcf ~{sample_id}.manta.vcf.gz
```

## Lumpy

### Cram to Bam conversion
`samtools view -h -T ~{ref_fasta} ~{input_cram} | samtools view -b -o ${inputBam}.bam`

### Get discordant reads
```sh
samtools view -h -@ ${threads} -F 1294 -u -b -h  ${inputBam} > ${inputBam}.temp && \
samtools sort -@ ${threads} -m 16G ${inputBam}.temp > discords.bam && rm ${inputBam}.temp
```

### Get split reads
```sh
samtools view -h -@ ${threads} ${bamToSplits} | \
/app/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
samtools view -@ ${threads} -b -u - > ${bamToSplits}.temp && \
samtools sort -@ ${threads} -m 16G ${bamToSplits}.temp > splits.bam && rm ${bamToSplits}.temp
```

### Run LumpyExpress
`lumpyexpress -B ${inputBam} -t ${threads} -S ${bamSplits} -D ${bamDiscords} -o ~{sample_name}_calls.vcf`


## PopIns2

### Link the reference genomes
`ln -s CHM13.fa genome.fa`
`ln -s CHM13.fa.fai genome.fa.fai`

### Assemble samples
```sh
popins2 assemble --sample sample1 sample1_CHM13.sorted.bam -t 24
...
popins2 assemble --sample sample10 sample10_CHM13.sorted.bam -t 24
```

### Merge
`popins2 merge -r PopIns2/CHM13/ -di`

### Contigmap
```sh
popins2 contigmap sample1 -t 24
...
popins2 contigmap sample10 -t 24
```

### Place
```sh
popins2 place-refalign
popins2 place-splitalign sample1
...
popins2 place-splitalign sample10
popins2 place-finish
```

### Genotype
```sh
popins2 genotype sample1
...
popins2 genotype sample10
```
