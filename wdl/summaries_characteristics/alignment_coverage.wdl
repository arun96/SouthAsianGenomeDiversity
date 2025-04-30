version 1.0
##
#WORKFLOW DEFINITION
workflow ComputeCoverageFlow {
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
  call ComputeCoverageTask{
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
    File output_coverage = ComputeCoverageTask.coverage
    File output_coverage_calc = ComputeCoverageTask.coverage_calc
  }

}

#Task Definitions
task ComputeCoverageTask {
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


  # Calls samtools to do the coverage calculations
  command <<<
   
    # Full summary file
    # samtools depth -a ~{input_cram} --reference ~{ref_fasta} > ~{sample_name}.coverage 
    
    # Samtools coverage - sort then call
    samtools sort ~{input_cram} --reference ~{ref_fasta} -o ~{sample_name}_sorted.cram --threads ~{num_threads}
    samtools coverage ~{sample_name}_sorted.cram --reference ~{ref_fasta} > ~{sample_name}.coverage
    

    # Calculated coverage
    samtools depth -a ~{input_cram} --reference ~{ref_fasta} | awk '{sum+=$3} END { print "Average = ",sum/NR}' > ~{sample_name}.calc.coverage
    
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
    File coverage = "~{sample_name}.coverage"
    File coverage_calc = "~{sample_name}.calc.coverage"
  }
}