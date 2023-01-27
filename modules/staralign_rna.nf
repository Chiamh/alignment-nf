// Using STAR aligner to align microbial reads to pangenome. Unlike bowtie2, reads are mapped as paired end because eukaryotic reads are not poly-cistronic.

process STARALIGN_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/STAR_out", mode: 'copy'
	
	input:
	path star_index
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_star_microbe_pangenome_aligned.bam"), emit: aligned
	tuple val(sample_id), path("${sample_id}_star_microbe_pangenome_aligned_filtered_cov.tsv"), emit: coverage
	tuple val(sample_id), path("*.out"), emit: logs
	tuple val(sample_id), path("${sample_id}_star_microbe_pangenome_unaligned_{1,2}.fastq.gz"), emit: unaligned
	
	when:
	params.map && params.process_rna
	
	script:

	"""
		
	STAR --runMode alignReads \\
		--runThreadN $task.cpus \\
		--outSAMtype BAM Unsorted \\
		--readFilesCommand zcat \\
		--genomeDir ${star_index} \\
		--outFileNamePrefix ${sample_id}. \\
		--readFilesIn ${reads_file[0]} ${reads_file[1]} \\
		--outReadsUnmapped Fastx
		
	if [ -f ${sample_id}.Unmapped.out.mate1 ]; then
	mv ${sample_id}.Unmapped.out.mate1 ${sample_id}_star_microbe_pangenome_unaligned_1.fastq
	gzip ${sample_id}_star_microbe_pangenome_unaligned_1.fastq
	fi
	if [ -f ${sample_id}.Unmapped.out.mate2 ]; then
	mv ${sample_id}.Unmapped.out.mate2 ${sample_id}_star_microbe_pangenome_unaligned_2.fastq
	gzip ${sample_id}_star_microbe_pangenome_unaligned_2.fastq
	fi
	
	if [ -f ${sample_id}.Aligned.out.bam ]; then
	mv ${sample_id}.Aligned.out.bam ${sample_id}_star_microbe_pangenome_aligned.bam
	fi
	
	star_helper.sh "${sample_id}"
			
	"""
}

