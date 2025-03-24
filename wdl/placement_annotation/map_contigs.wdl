version 1.0
##
#WORKFLOW DEFINITION
workflow MapContigsFlow {
  input {

    File reference_fasta
    File assembly_contigs

    String sample_name
    String output_ext = "1Kb_CHM13"

    String docker_image = "szarate/t2t_variants"
    Int preemptible_tries = 5
    Int threads = 16
    Int disk_size = 50
    Int machine_mem_size = 32
  }

  # assembles unaligned and poorly aligned reads 
  call MapContigsTask{
    input:
      reference_fasta = reference_fasta,
      output_ext = output_ext,
      assembly_contigs = assembly_contigs,
      sample_name = sample_name,
      docker_image = docker_image,
      preemptible_tries = preemptible_tries,
      num_threads = threads,
      disk_size = disk_size,
      machine_mem_size = machine_mem_size
  }

  #Outputs extracted reads and relevant files
  output {
    File alignment_file = MapContigsTask.alignment_sam
  }

}

#Task Definitions
task MapContigsTask {
  input {
    # Command parameters
    File reference_fasta
    File assembly_contigs

    String sample_name
    String output_ext

    # Runtime parameters
    String docker_image
    Int preemptible_tries
    Int num_threads

    Int disk_size
    Int machine_mem_size
  }

  Int runtime_memory = (machine_mem_size - 2)

  # Calls megahit to do the assembly
  command <<<

    minimap2 -a -I ~{runtime_memory}g -t ~{num_threads} ~{reference_fasta} ~{assembly_contigs} > ~{sample_name}_~{output_ext}_mapped.sam
    
  >>>

  #Run time attributes:
  runtime {
    docker: docker_image
    memory: machine_mem_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
    preemptible: preemptible_tries
    cpu: num_threads
  }
    
  #Outputs a BAM and BAI with the same sample name
  output {
    File alignment_sam = "${sample_name}_${output_ext}_mapped.sam"
  }
}