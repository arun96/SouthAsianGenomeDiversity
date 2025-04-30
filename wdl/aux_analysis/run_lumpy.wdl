version 1.0
##
#Task Definitions
task CramToBamTask {
	input {
		# Command parameters
		File ref_fasta
		File ref_fasta_index
		File ref_dict
		File input_cram
		String sample_name

		# Runtime parameters
		Int addtional_disk_size = 128 
		Int machine_mem_size = 32
		String docker_image
		Int preemptible_tries = 5
	}
		# Disk size
		Float output_bam_size = size(input_cram, "GB") / 0.60
		Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
		Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + addtional_disk_size


	#Calls samtools view to do the conversion
	command {
		set -eo pipefail
		samtools view -h -T ~{ref_fasta} ~{input_cram} | samtools view -b -o ~{sample_name}.bam
	}

	#Run time attributes:
	runtime {
		docker: docker_image
		memory: machine_mem_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
		
	#Outputs a BAM with the same sample name
	output {
		File outputBam = "~{sample_name}.bam"
	}
}

task getDiscordants{
	input {
		File inputBam
		Int threads
		String sample_name
		String docker_image

		# Runtime parameters
		Int disk_size = 128 
		Int mem_size = 32
		Int preemptible_tries = 5
	}

	command{
		samtools view -h -@ ${threads} -F 1294 -u -b -h  ${inputBam} > ${inputBam}.temp && \
		samtools sort -@ ${threads} -m 16G ${inputBam}.temp > discords.bam && rm ${inputBam}.temp
	}
	runtime{
		docker : docker_image
		memory: mem_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}

	output {
		File discordsBam="discords.bam"
	}
}

task getSplits{
	input {
		File bamToSplits
		Int threads
		String sample_name
		String docker_image

		# Runtime parameters
		Int disk_size = 128 
		Int mem_size = 32
		Int preemptible_tries = 5
	}

	command {
		samtools view -h -@ ${threads} ${bamToSplits} | \
		/app/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin | \
		samtools view -@ ${threads} -b -u - > ${bamToSplits}.temp && \
		samtools sort -@ ${threads} -m 16G ${bamToSplits}.temp > splits.bam && rm ${bamToSplits}.temp
	}
	runtime{
		docker : docker_image
		memory: mem_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output {
		File splitsBam="splits.bam"
	}
}

task lumpyexpress{
	input {
		File inputBam
		File bamSplits
		File bamDiscords
		Int threads
		String sample_name
		String docker_image

		# Runtime parameters
		Int disk_size = 128 
		Int mem_size = 32
		Int preemptible_tries = 5
	}

	command {
		lumpyexpress -B ${inputBam} -t ${threads} -S ${bamSplits} -D ${bamDiscords} -o ~{sample_name}_calls.vcf
	}
	runtime {
		docker : docker_image
		memory: mem_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: preemptible_tries
	}
	output{
		File outVCF="~{sample_name}_calls.vcf"
	}
}

workflow lumpyexpressFULL {
	input {
		Int threads
		String sample_name
		String docker_image = "erictdawson/lumpy-sv"

		# File inputBam
		File ref_fasta
		File ref_fasta_index
		File ref_dict
		File input_cram

		Int disk_size = 128
		Int mem_size = 32
		Int preemptible_tries = 5
	}

	call CramToBamTask{
		input:
			ref_fasta = ref_fasta,
			ref_fasta_index = ref_fasta_index,
			ref_dict = ref_dict,
			input_cram = input_cram,
			sample_name = sample_name,
			docker_image = docker_image
	}
   
	call getDiscordants{
		input:
			# inputBam=inputBam,
			inputBam = CramToBamTask.outputBam,
			threads=threads,
			sample_name=sample_name,
			docker_image=docker_image,
			disk_size = disk_size,
			mem_size = mem_size,
			preemptible_tries = preemptible_tries

	}

	call getSplits{
		input:
			# bamToSplits=inputBam,
			bamToSplits = CramToBamTask.outputBam,
			threads=threads,
			sample_name=sample_name,
			docker_image=docker_image,
			disk_size = disk_size,
			mem_size = mem_size,
			preemptible_tries = preemptible_tries

	}

	call lumpyexpress{
		input:
			inputBam=CramToBamTask.outputBam,
			threads=threads,
			sample_name=sample_name,
			docker_image=docker_image,
			bamSplits=getSplits.splitsBam,
			bamDiscords=getDiscordants.discordsBam,
			disk_size = disk_size,
			mem_size = mem_size,
			preemptible_tries = preemptible_tries
	}

	output {
		File lumpy_vcf = lumpyexpress.outVCF
	}
}