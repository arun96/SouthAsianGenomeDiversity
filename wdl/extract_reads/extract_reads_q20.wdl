version 1.0
##
#WORKFLOW DEFINITION
workflow ExtractReadsFlow {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name
    String gotc_docker = "szarate/t2t_variants"
    Int preemptible_tries = 5
    Int threads = 16
  }

  # pulls unaligned and poorly aligned reads 
  call ExtractReadsTask{
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      input_cram = input_cram,
      sample_name = sample_name,
      docker_image = gotc_docker,
      preemptible_tries = preemptible_tries,
      num_threads = threads
  }

  #Outputs extracted reads and relevant files
  output {
    # File output_q20 = ExtractReadsTask.q20
    File output_q20_sorted = ExtractReadsTask.q20_sorted
    File output_q20_R1 = ExtractReadsTask.q20_R1
    File output_q20_R2 = ExtractReadsTask.q20_R2
  }

}

#Task Definitions
task ExtractReadsTask {
  input {
    # Command parameters
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name

    # Runtime parameters
    Int addtional_disk_size = 50 
    Int machine_mem_size = 15
    String docker_image
    Int preemptible_tries
    Int num_threads
  }

    # TODO - fix the size of the VM
    Float output_size = size(input_cram, "GB") / 0.1
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_size + ref_size) + addtional_disk_size


  # Calls samtools to do the extractions
  command <<<

    samtools view -h -b --reference ~{ref_fasta} ~{input_cram} -q 20 -o ~{sample_name}_greater_than_q20.bam -U ~{sample_name}_q20.bam --threads ~{num_threads}
    samtools sort ~{sample_name}_q20.bam -o ~{sample_name}_q20_sorted.bam --threads ~{num_threads}
    samtools fastq ~{sample_name}_q20_sorted.bam -1 ~{sample_name}_q20_R1.fq -2 ~{sample_name}_q20_R2.fq --threads ~{num_threads}
    
  >>>

  #Run time attributes:
  #Use a docker with samtools. Set this up as a workspace attribute.
  #disk_size should equal input size + output size + buffer
  runtime {
    docker: docker_image
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    cpu: num_threads
  }
    
  #Outputs a BAM and BAI with the same sample name
  output {
    # File q20 = "~{sample_name}_q20.bam"
    File q20_sorted = "~{sample_name}_q20_sorted.bam"
    File q20_R1 = "~{sample_name}_q20_R1.fq"
    File q20_R2 = "~{sample_name}_q20_R2.fq"
  }
}