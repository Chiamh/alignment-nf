// Preprocessing AND decontamination or removal of human reads from metagenomes

process FASTP_DECONT_DNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/reads/DNA", mode: 'copy'
	
	input:
	path bwaidx_path
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") , emit: logs
	
	when:
	params.preprocess && params.decont && params.process_dna
	
	script:
	"""
	fastp -i ${reads_file[0]} -I ${reads_file[1]} --stdout -j ${sample_id}.json -h ${sample_id}.html | \\
	bwa mem -p -t $task.cpus ${bwaidx_path}/${params.bwaidx} - | \\
	samtools fastq -f12 -F256 -1 ${sample_id}_host_unmapped_1.fastq.gz -2 ${sample_id}_host_unmapped_2.fastq.gz -
		
	"""
}

