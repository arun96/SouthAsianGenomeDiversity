# Supplementary Guide 2: Assembly and Filtering

## Assembly of extracted reads using MEGAHIT:
`megahit --num-cpu-threads 48 --out-dir ./megahit -1 unmapped_R1.fq -2 unmapped_R2.fq -r mateMapped_R1.fq, mateMapped_R2.fq`

## Extracting only contigs >1 Kb:

### Turn multiline into single line per sequence
`awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $ASM > collapsed.fa`

### Extract 1Kb contigs
`awk -v RS='>[^\n]+\n' 'length() >= 1000 {printf "%s", prt $0} {prt = RT}' collapsed.fa > contigs_1Kb.fa`

## Running Centrifuge on the assembled contigs:
`centrifuge -x centrifuge-1.0.3-beta/DB/p_compressed+h+v --report-file /centrifuge.report -k 1 --host-taxids 9606 -f contigs_1Kb.fa > centrifuge.output`

## Running BLAST on the assembled contigs:
`blastn -db ref_prok_rep_genomes -query contigs_1Kb.fa -out blast_output.out`

The output files from Centrifuge and BLAST are then parsed to remove non-human sequence.