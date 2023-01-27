// Pangenome/gene catalog nucleotide alignment for metatranscriptomes (MTX)
// Reads are not mapped as paired end due to bacterial operons
/*
Requires a pre-built pangenome with a bowtie2 index
This process will concatenate the metagenomic R1 and R2 files before bowtie2 mapping to pangenome.
*/

process BT2ALIGN_RNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/bt2_out/RNA", mode: 'copy'
	
	input:
	path bt2_idx_path
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_aligned.bam"), emit: aligned
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_aligned_filtered_cov.tsv"), emit: coverage
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_unaligned.fastq.gz"), emit: unaligned
	tuple val(sample_id), path("${sample_id}_bt2.log"), emit: logs
	
	when:
	params.map && params.process_rna
	
	script:
	"""
	zcat ${reads_file[0]} ${reads_file[1]} | \\
	(bowtie2 -q -x ${bt2_idx_path}/${params.bt2_idx_name} -U - --un-gz "${sample_id}_bt2_microbe_pangenome_unaligned.fastq.gz" \\
	-p $task.cpus --very-sensitive) 2>"${sample_id}_bt2.log" | \\
	samtools view -bS - > "${sample_id}_bt2_microbe_pangenome_aligned.bam"

	bt2_helper.sh "${sample_id}"
	
	"""
}
