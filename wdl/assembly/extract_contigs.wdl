version 1.0
##
#WORKFLOW DEFINITION
workflow ExtractContigsFlow {
  meta {
        description: "Generate a summary file, and extract 1Kb contigs into another file."
  }
  
  input {
    File assembly_fasta
    String sample_name
    Int contig_size_kb = 1000
    String gotc_docker = "szarate/t2t_variants"
    Int preemptible_tries = 5
    Int threads = 1
  }

  # pulls unaligned and poorly aligned reads 
  call ExtractContigsTask{
    input:
      assembly_fasta = assembly_fasta,
      sample_name = sample_name,
      contig_size_kb = contig_size_kb,
      docker_image = gotc_docker,
      preemptible_tries = preemptible_tries,
      num_threads = threads
  }

  #Outputs extracted reads and relevant files
  output {
    # File output_collapsed_fasta = ExtractContigsTask.collapsed_fasta
    # File output_summary_file = ExtractContigsTask.summary_file
    File output_contigs_1kb_fasta = ExtractContigsTask.contigs_1kb_fasta
  }

}

#Task Definitions
task ExtractContigsTask {
  input {

    # Command Params
    File assembly_fasta
    String sample_name
    Int contig_size_kb

    # Runtime Params
    String docker_image
    Int preemptible_tries
    Int num_threads
    Int addtional_disk_size = 0
    Int machine_mem_size = 5

  }
    Int output_size = ceil(size(assembly_fasta, "GB") / 0.5)
    Int contig_size = contig_size_kb * 1000


  # Calls samtools to do the extractions
  command <<<

    # count contig size
    # cat  ~{assembly_fasta} | grep -c '.\{1000\}' > ~{sample_name}_contig_counts.txt
    # cat  ~{assembly_fasta} | grep -c '.\{5000\}' >> ~{sample_name}_contig_counts.txt
    # cat  ~{assembly_fasta} | grep -c '.\{10000\}' >> ~{sample_name}_contig_counts.txt

    # Extract 1Kb contigs
    # awk -v RS='>[^\n]+\n' 'length() >= 1000 {printf "%s", prt $0} {prt = RT}'  ~{assembly_fasta} > ~{sample_name}_contigs_1kb.fa
    awk -v n=~{contig_size} '/^>/{ if(l>n) print b; b=$0;l=0;next } {l+=length;b=b ORS $0}END{if(l>n) print b }' ~{assembly_fasta} > ~{sample_name}_contigs_~{contig_size_kb}kb.fa
    
  >>>

  #Run time attributes:
  #Use a docker with samtools. Set this up as a workspace attribute.
  #disk_size should equal input size + output size + buffer
  runtime {
    docker: docker_image
    memory: machine_mem_size + " GB"
    disks: "local-disk " + output_size + " HDD"
    preemptible: preemptible_tries
    cpu: num_threads
  }
    
  #Outputs a BAM and BAI with the same sample name
  output {
    # File collapsed_fasta = "~{sample_name}_collapsed.fa"
    # File summary_file = "~{sample_name}_contig_counts.txt"
    File contigs_1kb_fasta = "~{sample_name}_contigs_~{contig_size_kb}kb.fa" 
  }
}