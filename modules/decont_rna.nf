// decontamination/removal of human reads using STAR aligner to get microbial reads
// To mix nextflow and shell variables: https://www.nextflow.io/docs/latest/process.html#process-shell

process DECONT_RNA {
	label "process_high"
	tag "${sample_id}"
	if(params.save_intermediates && params.remove_rRNA || params.save_intermediates && params.dedupe){
		publishDir "${params.outdir}/reads/RNA/temp", mode: 'copy', pattern: '*.fastq.gz'
	} else if (!params.remove_rRNA && !params.dedupe){
		publishDir "${params.outdir}/reads/RNA",mode: 'copy', pattern: '*.fastq.gz'
	}
	
	input:
	path decont_star_index
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_host_unmapped_{1,2}.fastq.gz"), emit: reads
	
	when:
	!params.preprocess && params.decont && params.process_rna
	
	script:

	"""
		
	STAR --runMode alignReads \\
		--runThreadN $task.cpus \\
		--outSAMtype None \\
		--readFilesCommand zcat \\
		--genomeDir ${decont_star_index} \\
		--outFileNamePrefix ${sample_id}. \\
		--readFilesIn ${reads_file[0]} ${reads_file[1]} \\
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

