version 1.0

workflow runMinigraphMap {

    input {
        File reads1
        File reads2
        File pangenome_graph

        String sample_name
        String args
        String output_file_tag = "minigraph"
        String docker_image = "humanpangenomics/hpp_minigraph:latest"
        Int preemptible_tries = 5
        Int num_threads = 16
        Int disk_size = 64
        Int machine_mem_size = 32
    }
    call minigraphMap {
        input:
            reads1 = reads1,
            reads2 = reads2,
            pangenome_graph = pangenome_graph,

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
        File output_gaf_gz = minigraphMap.output_gaf_gz
    }
}



task minigraphMap {

    input {
        File reads1
        File reads2
        File pangenome_graph

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
    String output_gaf_gz = "${sample_name}.${output_file_tag}.gaf.gz"

    command <<<

        # set -eux -o pipefail

        reads1_unzipped=$(basename -- "~{reads1}")
        ## first check if reads1 needs to be unzipped
        if [[ $reads1_unzipped =~ \.gz$ ]]; then
            cp ~{reads1} .
            gunzip -f $reads1_unzipped
            reads1_unzipped="${reads1_unzipped%.gz}"
        else
            ln -s ~{reads1}
        fi

        reads2_unzipped=$(basename -- "~{reads2}")
        ## first check if reads2 needs to be unzipped
        if [[ $reads2_unzipped =~ \.gz$ ]]; then
            cp ~{reads2} .
            gunzip -f $reads2_unzipped
            reads2_unzipped="${reads2_unzipped%.gz}"
        else
            ln -s ~{reads2}
        fi

        # zipped version
        minigraph -t ~{num_threads} ~{args} ~{pangenome_graph} ${reads1_unzipped} ${reads2_unzipped} | gzip > ~{output_gaf_gz}

    >>>

    output {

        File output_gaf_gz  = output_gaf_gz
    }

    runtime {
        memory: machine_mem_size + " GB"
        cpu: num_threads
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible_tries
    }
}