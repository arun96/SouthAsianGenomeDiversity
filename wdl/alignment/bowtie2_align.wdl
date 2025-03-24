version 1.0

workflow bowtie2_align {

    input {
        File reads1
        File reads2
        File reference_fasta

        String sample_name
        String args = ""
        String output_file_tag = "CHM13v2"
        String docker_image = "mschatz/wga-essentials"
        Int preemptible_tries = 5
        Int num_threads = 32
        Int disk_size = 1000
        Int machine_mem_size = 64
    }
    call bowtie2_align_task {
        input:
            reads1 = reads1,
            reads2 = reads2,
            reference_fasta = reference_fasta,

            sample_name = sample_name,
            args = args,
            output_file_tag = output_file_tag,

            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            num_threads = num_threads,
            disk_size = disk_size,
            machine_mem_size = machine_mem_size
    }

    output {
        File alignment_bam = bowtie2_align_task.alignment_bam
    }
}


task bowtie2_align_task {

    input {
        File reads1
        File reads2
        File reference_fasta

        String sample_name
        String args
        String output_file_tag

        String docker_image
        Int preemptible_tries
        Int num_threads
        Int disk_size
        Int machine_mem_size
    }

    # zipped version
    String index_files = "./${output_file_tag}"

    command <<<

       bowtie2-build ~{reference_fasta} ~{index_files} -p ~{num_threads}

       bowtie2 -p ~{num_threads} -x ~{index_files} -1 ~{reads1} -2 ~{reads2} | samtools view -bS - > ~{sample_name}_~{output_file_tag}.bam

    >>>

    output {

        File alignment_bam  = "~{sample_name}_~{output_file_tag}.bam"
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}