version 1.0
##
#WORKFLOW DEFINITION
workflow ExtractMapQFlow {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_cram
    String sample_name
    String gotc_docker = "quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    Int preemptible_tries = 5
    Int threads = 16
  }

  # Computes coverage
  call ExtractMapQTask{
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

  # Outputs a processed version
  output {
    File output_mapq_list = ExtractMapQTask.mapq_list
  }

}

#Task Definitions
task ExtractMapQTask {
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


  # Calls samtools
  command <<<

    # View and extract mapq field
    # samtools view ~{input_cram} --reference ~{ref_fasta} | cut -f5 > ~{sample_name}.mapq.txt
    # samtools view ~{input_cram} --reference ~{ref_fasta} | cut -f5 > ~{sample_name}.tmp.mapq.txt
    # cat ~{sample_name}.tmp.mapq.txt | sort | uniq -c | sort -k2n > ~{sample_name}.mapq.txt
    
    samtools view ~{input_cram} --reference ~{ref_fasta} | cut -f5 > ~{sample_name}.mapq.txt
    
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
    
  # Output files
  output {
    File mapq_list = "~{sample_name}.mapq.txt"
  }
}