version 1.0

workflow bedtools_intersect_genome {

    input {
        File alignment_bed
        File annotation_file

        String sample_name

        String output_file_tag = "contigs_unaligned_1kb_CHM13_mapped"

        String docker_image = "mschatz/wga-essentials"

        Int preemptible_tries = 5
        Int num_threads = 32
        Int disk_size = 1000
        Int machine_mem_size = 64
    }
    call bedtools_intersect_genome_task {
        input:
            alignment_bed = alignment_bed,
            annotation_file = annotation_file,

            sample_name = sample_name,

            output_file_tag = output_file_tag,

            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            num_threads = num_threads,
            disk_size = disk_size,
            machine_mem_size = machine_mem_size
    }

    output {
        File sorted_bed_file = bedtools_intersect_genome_task.sorted_bed_file
    }
}


task bedtools_intersect_genome_task {

    input {
        File alignment_bed
        File annotation_file

        String sample_name

        String output_file_tag

        String docker_image
        Int preemptible_tries
        Int num_threads
        Int disk_size
        Int machine_mem_size
    }

    command <<<

       bedtools intersect -a ~{alignment_bed} -b ~{annotation_file} -wb > ~{sample_name}_~{output_file_tag}.intersect.bed

    >>>

    output {

        File sorted_bed_file  = "~{sample_name}_~{output_file_tag}.intersect.bed"
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}