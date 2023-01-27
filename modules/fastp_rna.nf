// Preprocessing reads with fastp. 

process FASTP_RNA {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/reads/RNA", mode: 'copy', pattern: '*.{html,json}'
	if(params.save_intermediates && params.remove_rRNA || params.save_intermediates && params.dedupe){
		publishDir "${params.outdir}/reads/RNA/temp", mode: 'copy', pattern: '*.fastq.gz'
	} else if (!params.remove_rRNA && !params.dedupe){
		publishDir "${params.outdir}/reads/RNA",mode: 'copy', pattern: '*.fastq.gz'
	}
	
	
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

