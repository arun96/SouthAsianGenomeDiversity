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
    File output_unmapped_R1 = ExtractReadsTask.unmapped_R1
    File output_unmapped_R2 = ExtractReadsTask.unmapped_R2
    File output_mateMapped_R1 = ExtractReadsTask.mateMapped_R1
    File output_mateMapped_R2 = ExtractReadsTask.mateMapped_R2
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
    Int addtional_disk_size = 20 
    Int machine_mem_size = 15
    String docker_image
    Int preemptible_tries
    Int num_threads
  }
    Float output_size = size(input_cram, "GB") / 0.50
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_cram, "GB") + output_size + ref_size) + addtional_disk_size


  # Calls samtools to do the extractions
  command <<<
   
    samtools fastq -f 12 --reference ~{ref_fasta} ~{input_cram} -1 ~{sample_name}_unmapped_R1.fq -2 ~{sample_name}_unmapped_R2.fq --threads ~{num_threads}

    samtools fastq -f 68 -F 8 --reference ~{ref_fasta} ~{input_cram} --threads ~{num_threads} > ~{sample_name}_mateMapped_R1.fq

    samtools fastq -f 132 -F 8 --reference ~{ref_fasta} ~{input_cram} --threads ~{num_threads} > ~{sample_name}_mateMapped_R2.fq
    
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
    File unmapped_R1 = "~{sample_name}_unmapped_R1.fq"
    File unmapped_R2 = "~{sample_name}_unmapped_R2.fq"
    File mateMapped_R1 = "~{sample_name}_mateMapped_R1.fq"
    File mateMapped_R2 = "~{sample_name}_mateMapped_R2.fq"
  }
}