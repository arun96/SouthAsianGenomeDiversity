version 1.0
# Inspired by: https://dockstore.org/workflows/github.com/human-pangenomics/hpp_production_workflows/MapBlocks:asset?tab=info 

workflow generate_alignment_bed {

    input {
        File alignment_file

        String sample_name

        String output_file_tag = "contigs_unaligned_1kb_CHM13_mapped"

        String docker_image = "mschatz/wga-essentials"

        Int preemptible_tries = 5
        Int num_threads = 32
        Int disk_size = 1000
        Int machine_mem_size = 64
    }
    call generate_alignment_bed_task {
        input:
            alignment_file = alignment_file,

            sample_name = sample_name,

            output_file_tag = output_file_tag,

            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            num_threads = num_threads,
            disk_size = disk_size,
            machine_mem_size = machine_mem_size
    }

    output {
        File sorted_bed_file = generate_alignment_bed_task.sorted_bed_file
    }
}


task generate_alignment_bed_task {

    input {
        File alignment_file

        String sample_name

        String output_file_tag

        String docker_image
        Int preemptible_tries
        Int num_threads
        Int disk_size
        Int machine_mem_size
    }

    command <<<

       samtools view -Sb ~{alignment_file} | bedtools bamtobed -i - | bedtools sort -i - > ~{sample_name}_~{output_file_tag}.sorted.bed

    >>>

    output {

        File sorted_bed_file  = "~{sample_name}_~{output_file_tag}.sorted.bed"
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}