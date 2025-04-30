version 1.0

workflow jellyfish_count {

    input {
        File reads1
        File reads2

        String sample_name
        
        String args = ""
        
        String docker_image = "mschatz/wga-essentials"
        Int preemptible_tries = 5
        Int num_threads = 32
        Int disk_size = 1000
        Int machine_mem_size = 64
    }
    call jellyfish_count_task {
        input:
            reads1 = reads1,
            reads2 = reads2,

            sample_name = sample_name,
            
            args = args,

            docker_image = docker_image,
            preemptible_tries = preemptible_tries,
            num_threads = num_threads,
            disk_size = disk_size,
            machine_mem_size = machine_mem_size
    }

    output {
        File mer_counts = jellyfish_count_task.mer_counts
        File mer_histo = jellyfish_count_task.mer_histo
    }
}


task jellyfish_count_task {

    input {
        File reads1
        File reads2

        String sample_name
        
        String args

        String docker_image
        Int preemptible_tries
        Int num_threads
        Int disk_size
        Int machine_mem_size
    }

    command <<<
       jellyfish count -C -m 21 -s 100G -t 12 <(zcat ~{reads1}) <(zcat ~{reads2}) -o ~{sample_name}_mer_counts.jf

       jellyfish histo ~{args} -t 12 ~{sample_name}_mer_counts.jf > ~{sample_name}_mer_counts.histo

    >>>

    output {

        File mer_counts  = "~{sample_name}_mer_counts.jf"
        File mer_histo = "~{sample_name}_mer_counts.histo"
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}