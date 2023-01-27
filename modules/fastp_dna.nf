// Preprocessing reads with fastp. 

process FASTP_DNA {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/reads/DNA", mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_fastp_{1,2}.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") , emit: logs
	
	when:
	params.preprocess && !params.decont
	
	script:
	"""
	fastp -i ${reads_file[0]} -I ${reads_file[1]} \\
		--out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \\
		-j ${sample_id}.json -h ${sample_id}.html		
	"""
}

