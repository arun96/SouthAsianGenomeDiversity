version 1.0
##
#WORKFLOW DEFINITION
workflow MegahitUnalignedFlow {
  input {

    File unmapped_R1
    File unmapped_R2
    File mateMapped_R1
    File mateMapped_R2

    String sample_name

    String docker_image = "vout/megahit"
    Int preemptible_tries = 5
    Int threads = 16
    Int disk_size = 50
    Int machine_mem_size = 16
  }

  # assembles unaligned and poorly aligned reads 
  call MegahitUnalignedTask{
    input:
      unmapped_R1 = unmapped_R1,
      unmapped_R2 = unmapped_R2,
      mateMapped_R1 = mateMapped_R1,
      mateMapped_R2 = mateMapped_R2,
      sample_name = sample_name,
      docker_image = docker_image,
      preemptible_tries = preemptible_tries,
      num_threads = threads,
      disk_size = disk_size,
      machine_mem_size = machine_mem_size
  }

  #Outputs extracted reads and relevant files
  output {
    File contigs = MegahitUnalignedTask.final_contigs
    File log = MegahitUnalignedTask.assembly_log
  }

}

#Task Definitions
task MegahitUnalignedTask {
  input {
    # Command parameters
    File unmapped_R1
    File unmapped_R2
    File mateMapped_R1
    File mateMapped_R2

    String sample_name

    # Runtime parameters
    String docker_image
    Int preemptible_tries
    Int num_threads

    Int disk_size
    Int machine_mem_size
  }

  Int assembly_memory = (machine_mem_size - 1) * 1000000000

  # Calls megahit to do the assembly
  command <<<

    megahit -1 ~{unmapped_R1} -2 ~{unmapped_R2} -r ~{mateMapped_R1},~{mateMapped_R2} -t ~{num_threads} -m ~{assembly_memory} -o assembly

    mv assembly/final.contigs.fa assembly/~{sample_name}_final.contigs.fa
    
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
    File final_contigs = "assembly/${sample_name}_final.contigs.fa"
    File assembly_log = "assembly/log"
  }
}