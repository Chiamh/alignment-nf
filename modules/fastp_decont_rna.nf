// Preprocessing AND decontamination/removal of human reads using STAR aligner to get microbial reads
// To mix nextflow and shell variables: https://www.nextflow.io/docs/latest/process.html#process-shell

process FASTP_DECONT_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/reads/RNA", mode: 'copy', pattern: '*.{html,json}'
	if(params.save_intermediates && params.remove_rRNA || params.save_intermediates && params.dedupe){
		publishDir "${params.outdir}/reads/RNA/temp", mode: 'copy', pattern: '*.fastq.gz'
	} else if (!params.remove_rRNA && !params.dedupe){
		publishDir "${params.outdir}/reads/RNA",mode: 'copy', pattern: '*_unmapped_{1,2}.fastq.gz'
	}
	
	input:
	path decont_star_index
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_host_unmapped_{1,2}.fastq.gz"), emit: reads
	tuple path("${sample_id}.html"), path("${sample_id}.json") emit: logs
	tuple val(sample_id), path("${sample_id}_fastp_{1,2}.fastq.gz"), emit: fastpreads
	
	when:
	params.preprocess && params.decont && params.process_rna
	
	script:

	"""
	fastp -i ${reads_file[0]} -I ${reads_file[1]} \\
	--out1 ${sample_id}_fastp_1.fastq.gz --out2 ${sample_id}_fastp_2.fastq.gz \\
	-j ${sample_id}.json -h ${sample_id}.html
		
	STAR --runMode alignReads \\
		 --runThreadN $task.cpus \\
		 --outSAMtype None \\
		 --readFilesCommand zcat \\
		 --genomeDir ${decont_star_index} \\
		 --outFileNamePrefix ${sample_id}. \\
		 --readFilesIn ${sample_id}_fastp_1.fastq.gz ${sample_id}_fastp_2.fastq.gz \\
		 --outReadsUnmapped Fastx
		
	if [ -f ${sample_id}.Unmapped.out.mate1 ]; then
	mv ${sample_id}.Unmapped.out.mate1 ${sample_id}_host_unmapped_1.fastq
	gzip ${sample_id}_host_unmapped_1.fastq
	fi
	if [ -f ${sample_id}.Unmapped.out.mate2 ]; then
	mv ${sample_id}.Unmapped.out.mate2 ${sample_id}_host_unmapped_2.fastq
	gzip ${sample_id}_host_unmapped_2.fastq
	fi
	"""
}

