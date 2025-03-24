version 1.0

workflow fastqc_analysis {

    input {
        File reads_file
        String reads_extension = ".fastq.gz"

        String sample_name
        
        String docker_image = "mschatz/wga-essentials"
        Int preemptible_tries = 5
        Int num_threads = 16
        Int disk_size = 1000
        Int machine_mem_size = 64
    }
    call fastqc_analysis_task {
        input:
            reads_file = reads_file,
            reads_extension = reads_extension,

            sample_name = sample_name,

            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            num_threads = num_threads,
            disk_size = disk_size,
            machine_mem_size = machine_mem_size
    }

    output {
        File fastqc_html = fastqc_analysis_task.fastqc_html
        File fastqc_zip = fastqc_analysis_task.fastqc_zip
    }
}


task fastqc_analysis_task {

    input {
        File reads_file
        String reads_extension
        String sample_name

        String docker_image
        Int preemptible_tries
        Int num_threads
        Int disk_size
        Int machine_mem_size
    }
    String output_name = basename(reads_file, reads_extension)

    command <<<

       fastqc -t ~{num_threads} --noextract -o ./ ~{reads_file}

    >>>

    output {

        File fastqc_html  = "~{output_name}_fastqc.html"
        File fastqc_zip = "~{output_name}_fastqc.zip"
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}