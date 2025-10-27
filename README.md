# South Asian Genome Diversity

This is our repository containing all links and materials relevant to our work on investigating the diversity in and the representation of South Asians in genomic datasets. This repository will contain links to all talks, posters and papers relevant to this work, as well as the analysis scripts and tools used in this project.

This work was posted to bioRxiv in May 2025, and can be found here: ["Assembling unmapped reads reveals hidden variation in South Asian genomes", Das et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.05.14.653340v1.abstract).

We have recently completed the first round of revisions and the paper is once again under review, and all updated versions of the manuscript will be uploaded to bioRxiv.

## Repository Layout

All links relevant to this project (manuscript, past posters, workspaces and workflows, and data we generated) can found in the "Links" section of this README.

In the `wdl` folder, you can find all relevant WDL files used throughout our pipeline for aligning input reads, extracting unmapped reads, assembling and filtering the unmapped read contigs, placement against reference genomes, intersection with annotated regions and other `.bed` files, and any other miscellaneous scripts. These WDL files are divided into the relevant subfolders.

In the `scripts` folder, you can find any relevant Python notebooks and analysis scripts from local or cluster analysis.


### WDL Workflow Descriptions

`alignment/`
- `alignment/bowtie2_align.wdl`: Align a pair of input read files against a chosen reference fasta, and output a .bam alignment file.
- `alignment/minigraph_align.wdl`: Align a pair of input read files against a chosen minigraph pangenome graph, and output a .gaf.gz alignment file.

`extract_reads/`
- `extract_reads/extract_reads_unaligned.wdl`: Extracts four sets of unmapped reads from an alignment cram file, provided the reference fasta.
- `extract_reads/extract_reads_q20.wdl`: Extracts a pair of files containing poorly aligned reads from an alignment cram file, provided the reference fasta.
- `extract_reads/extract_reads_bam_unaligned.wdl`:  Extracts four sets of unmapped reads from an alignment bam file.
- `extract_reads/extract_reads_bam.wdl`: Extracts four sets of unmapped reads AND a pair of files containing poorly aligned reads from an alignment bam file.
- `extract_reads/extract_reads`: Extracts four sets of unmapped reads AND a pair of files containing poorly aligned reads from an alignment cram file, provided the reference fasta.

`assembly/`
- `assembly/extract_contigs.wdl`: Extracts all contigs above a chosen size into a separate .fa file.
- `assembly/megahit_unaligned.wdl`: Assemble all provided unmapped reads into contigs using default MEGAHIT parameters.
- `assembly/megahit_q20.wdl`: Assemble a pair of poorly aligned read files into contigs using default MEGAHIT parameters.

`aux_analysis/`
- `aux_analysis/fastqc_analysis.wdl`: Runs fastqc on reads files, generating the standard output.
- `aux_analysis/jellyfish_count.wdl`: Runs jellyfish on a pair of read files, outputing the counts and histogram/
- `aux_analysis/run_manta.wdl`: Runs manta on an alignment file - adapted from GATK-SV's Manta.wdl.
- `aux_analysis/run_lumpy.wdl`: Runs lumpy on an alignment file.

`placement_annotation/`
- `placement_annotation/map_contigs.wdl`: Maps input contigs to a chosen reference genome using Minimap2, outputting a sam alignment file.
- `placement_annotation/generate_alignment_bed.wdl`: Generates a bed file from an input contig sam format alignment file.
- `placement_annotation/bedtools_intersect_genome.wdl`: Intersects an alignment bed file against a provided annotation bed file, using the `-wb` parameter.

`summaries_characteristics/`
- `summaries_characteristics/alignment_coverage.wdl`: Output of the `samtools coverage` call fora chosen cram alignment file.
- `summaries_characteristics/alignment_coverage_bam.wdl`: Output of the `samtools coverage` call fora chosen bam alignment file.
- `summaries_characteristics/extract_contigs_summarize.wdl`: Generates an assembly summary file, and outputs all >1Kb contigs to a new fasta file.
- `summaries_characteristics/extract_fasta_summary.wdl`: Outputs seelcted sections of the results of `assembly-stats`.
- `summaries_characteristics/extract_mapq.wdl`: Outputs the mapping qualities from an input cram file.
- `summaries_characteristics/extract_megahit_summary.wdl`: Outputs the summary from a MEGAHIT assembly.
- `summaries_characteristics/extract_total_sequence_length.wdl`: Outputs the total sequence length of an input fasta file.


## Links

This work was posted to bioRxiv in May 2025, and can be found here: ["Assembling unmapped reads reveals hidden variation in South Asian genomes", Das et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.05.14.653340v1.abstract).

A selected subset of the key data from our work can be found on [Zenodo](https://zenodo.org/records/17419004). This repository contains the assembled contigs against CHM13, the results of collapsing sequence down based on sharedness/uniqueness, RepeatMasker analysis, and gene and GWAS intersections.

Our Anvil workspace can be found [here](https://anvil.terra.bio/#workspaces/anvil-dash-research/south-asian-genome). Most of our larger analyses were carried out on Anvil, and the scripts used can be found here, as well as the data generated by our work.

The supplementary tables referenced in the manuscript can be accessed [here](https://drive.google.com/drive/folders/1b2c1Xycqyr8_DPuFImsXIWtEAZlALJjy?usp=sharing). These tables include details of all genes intersected by the SAS and non-SAS placed contigs, and the top BLAST results of the unplaced contigs and the high scoring RNA-alignment contigs.

Full size versions of all figures used in the manuscript and supplement can also be found in the same [folder](https://drive.google.com/drive/folders/1b2c1Xycqyr8_DPuFImsXIWtEAZlALJjy?usp=sharing). The figures are split into folders depending on if they are in the main manuscript or the supplement, and lossless SVG versions should be available for all main figures and most supplemental ones.


## Past Presentations and Posters of this work

[Biology of Genomes 2025 Poster](https://drive.google.com/file/d/1a_KKooTs3Xxs78HY3Cl57SwL7QpgGZbP/view?usp=sharing) - shown at the Biology of Genomes 2025 meeting @ CSHL.

[BIODATA 2024 Poster](https://drive.google.com/file/d/1b5QrxCZGZ1VqBroGnVwDCJPVRb6DikzO/view?usp=sharing) - shown at the Biological Data Science 2024 Meeting @ CSHL

[GI 2023 Poster](https://drive.google.com/file/d/1sMByc7kWJl9ipy14WLtpWlJ2mlwt0hOc/view?usp=sharing) - shown at the Genome Informatics 2023 Meeting @ CSHL.

[ASHG 2023 Poster](https://drive.google.com/file/d/1N4YAV44Velab-i63iTrwN8bNHwAhO975/view?usp=sharing) - shown at the ASHG 2023 Meeting @ Washington DC.

## Contact

If you're coming here from my poster, and want to leave a comment/note for me, please do so using the following form: [Arun Das SA Genomes Contact Form](https://docs.google.com/forms/d/1rbrBBXTupyRSqqrBFbuvgqbC6mOgn53jzAnOVLSAfxY/)

Otherwise, please send me an email or open an issue on this repository.

Thank you!
