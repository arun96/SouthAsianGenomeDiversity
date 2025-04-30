version 1.0
##
#WORKFLOW DEFINITION
workflow ExtractReadsFlow {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_bam
    String sample_name
    String gotc_docker = "mschatz/wga-essentials"
    Int preemptible_tries = 5
    Int threads = 16
  }

  # pulls unaligned and poorly aligned reads 
  call ExtractReadsTask{
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      input_bam = input_bam,
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
    File input_bam
    String sample_name

    # Runtime parameters
    Int addtional_disk_size = 20 
    Int machine_mem_size = 15
    String docker_image
    Int preemptible_tries
    Int num_threads
  }
    Float output_bam_size = size(input_bam, "GB") / 0.10
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
    Int disk_size = ceil(size(input_bam, "GB") + output_bam_size + ref_size) + addtional_disk_size


  # Calls samtools to do the extractions
  command {
    set -eo pipefail

    samtools fastq -f 12 ~{input_bam} -1 ~{sample_name}_unmapped_R1.fq -2 ~{sample_name}_unmapped_R2.fq --threads ~{num_threads}

    samtools fastq -f 68 -F 8 ~{input_bam} --threads ~{num_threads} > ~{sample_name}_mateMapped_R1.fq 

    samtools fastq -f 132 -F 8 ~{input_bam} --threads ~{num_threads} > ~{sample_name}_mateMapped_R2.fq

    # samtools view -b -f 4 ~{input_bam} --threads ~{num_threads} > ~{sample_name}_unmapped.bam
    # samtools sort ~{sample_name}_unmapped.bam -o ~{sample_name}_unmapped_sorted.bam --threads ~{num_threads}

    # samtools view -h -b ~{input_bam} -q 20 -U ~{sample_name}_q20.bam --threads ~{num_threads}
    # samtools sort ~{sample_name}_q20.bam -o ~{sample_name}_q20_sorted.bam --threads ~{num_threads}
    # samtools fastq ~{sample_name}_q20_sorted.bam -1 ~{sample_name}_q20_R1.fq -2 ~{sample_name}_q20_R2.fq --threads ~{num_threads}
  }

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