// extract non-host RNA reads from bam files (alignment to host genome)

process BAM_TO_FASTQ_RNA {
	label "process_750GB"
	tag "${sample_id}"
	publishDir "${params.outdir}/raw/RNA", mode: 'copy'
	publishDir "${params.outdir}/raw/non_host_singletons/RNA", mode: 'copy', pattern: '_singleton.fastq.gz'
	
	input:
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple val(sample_id), path("${sample_id}_host_unmapped_singleton.fastq.gz"), emit: singletons
	tuple val(sample_id), path("${sample_id}_flagstat.txt"), emit: flagstat
	
	when:
	params.process_rna
	
	script:
	"""
	samtools flagstat -O tsv -@ $task.cpus ${reads_file[0]} > ${sample_id}_flagstat.tsv
	
	samtools collate -u -O --threads $task.cpus ${reads_file[0]} | \\
	samtools fastq -f12 -F256 -1 ${sample_id}_host_unmapped_1.fastq.gz -2 ${sample_id}_host_unmapped_2.fastq.gz -s ${sample_id}_host_unmapped_singleton.fastq.gz -
	
	"""
}

