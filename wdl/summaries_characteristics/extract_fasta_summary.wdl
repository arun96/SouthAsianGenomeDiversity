version 1.0

workflow extract_fasta_summary {
    meta {
        description: "Extract the assembly-stats summary of a fasta file"
    }

    input {
        File seq_file
        Int preemptible_tries = 3
        Int threads = 1
        String gotc_docker = "sangerpathogens/assembly-stats"
    }

    call extract_fasta_summary_line {
        input:
            seq_file      = seq_file,
            docker_image = gotc_docker,
            preemptible_tries = preemptible_tries,
            num_threads = threads
    }

    output {
        String summary_line = extract_fasta_summary_line.summary
    }
}

task extract_fasta_summary_line {
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
        assembly-stats ~{seq_file} | awk 'NR == 2 || NR==3' | tr '\n' ' '
    }
    

    runtime {
        docker: docker_image
    }

    output {
        String summary = read_string(stdout())
    }
}