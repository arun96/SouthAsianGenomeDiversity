version 1.0

workflow extract_total_sequence_length {
    meta {
        description: "Extract the total amount of sequence in a fasta file"
    }

    input {
        File seq_file
        Int preemptible_tries = 3
        Int threads = 1
        String gotc_docker = "mschatz/wga-essentials"
    }

    call extract_total_sequence_length_line {
        input:
            seq_file      = seq_file,
            docker_image = gotc_docker,
            preemptible_tries = preemptible_tries,
            num_threads = threads
    }

    output {
        String summary_line = extract_total_sequence_length_line.summary
    }
}

task extract_total_sequence_length_line {
    input {
        File seq_file
        
        # Runtime Params
        Int addtional_disk_size = 1 
        Int machine_mem_size = 1
        String docker_image
        Int preemptible_tries
        Int num_threads
    }

    command {
        grep -v ">" ~{seq_file} | tr -d '\n' | wc -c
    }
    

    runtime {
        docker: docker_image
    }

    output {
        String summary = read_string(stdout())
    }
}