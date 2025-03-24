version 1.0

workflow extract_megahit_summary {
    meta {
        description: "Extract the final summary of a megahit assembly run."
    }

    input {
        File log_file
        Int preemptible_tries = 3
        Int threads = 1
        String gotc_docker = "mschatz/wga-essentials"
    }

    call extract_megahit_summary_line {
        input:
            log_file      = log_file,
            docker_image = gotc_docker,
            preemptible_tries = preemptible_tries,
            num_threads = threads
    }

    output {
        String summary_line = extract_megahit_summary_line.summary
    }
}

task extract_megahit_summary_line {
    input {
        File log_file
        
        # Runtime Params
        Int addtional_disk_size = 1 
        Int machine_mem_size = 1
        String docker_image
        Int preemptible_tries
        Int num_threads
    }

    command {
        tail -2 ~{log_file} | head -1 | cut -c 23-
    }
    

    runtime {
        docker: docker_image
    }

    output {
        String summary = read_string(stdout())
    }
}