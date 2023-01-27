// rRNA removal (human, fungal and bacteria) from RNAseq data using bbduk.sh
// If params.dedupe = false, then rRNA removal is the "final" step of decontamination
//The publishDir can be included within a Groovy condition statement.

process RIBOFILTER {
	label "process_medium"
	tag "${sample_id}"
	
	if(params.save_intermediates && params.dedupe){
		publishDir "${params.outdir}/reads/RNA/temp", mode: 'copy', pattern: '*.fastq.gz'
	} else if (!params.dedupe){
		publishDir "${params.outdir}/reads/RNA",mode: 'copy', pattern: '*.fastq.gz'
	}
	publishDir "${params.outdir}/reads/RNA", mode: 'copy', pattern: '*.log'
	
	input:
	path ribokmers
	tuple val(sample_id), path(reads)
	
	output:
	tuple val(sample_id), path("${sample_id}*_{1,2}.fastq.gz"), emit: reads
	tuple val(sample_id), path("*_rRNAfilter.log") , emit: logs
	
	when:
	params.remove_rRNA && params.process_rna
	
	script:
	"""
	bbduk.sh in=${reads[0]} in2=${reads[1]} \\
	out=${sample_id}_mRNA_1.fastq.gz out2=${sample_id}_mRNA_2.fastq.gz \\
	k=31 \\
	ref=${ribokmers} stats=${sample_id}_rRNAfilter.log
	"""

	
	
