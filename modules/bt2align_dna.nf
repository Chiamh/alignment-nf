// Pangenome/gene catalog nucleotide alignment for metagenomes (MGX)

/*
Requires a pre-built pangenome with a bowtie2 index
This process will concatenate the metagenomic R1 and R2 files before bowtie2 mapping to pangenome by default unless force_paired == true
*/

process BT2ALIGN_DNA {
	label "process_high"
	tag "${sample_id}"
	publishDir "${params.outdir}/bt2_out/DNA", mode: 'copy'
	
	input:
	path bt2_idx_path
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_aligned.bam"), emit: aligned
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_aligned_filtered_cov.tsv"), emit: coverage
	tuple val(sample_id), path("${sample_id}_bt2_microbe_pangenome_unaligned.fastq.gz"), emit: unaligned
	tuple val(sample_id), path("${sample_id}_bt2.log"), emit: logs
	
	when:
	params.map && params.process_dna
	
	script:
	if (params.force_paired){

	"""
        (bowtie2 -q -x ${bt2_idx_path}/${params.bt2_idx_name} -1 ${reads_file[0]} -2 ${reads_file[1]} --un-gz "${sample_id}_bt2_microbe_pangenome_unaligned.fastq.gz" \\
        -p $task.cpus --very-sensitive) 2>"${sample_id}_bt2.log" | \\
        samtools view -bS - > "${sample_id}_bt2_microbe_pangenome_aligned.bam"

        bt2_helper_paired.sh "${sample_id}"

        """
	} else if (params.force_paired == false){
	"""
	zcat ${reads_file[0]} ${reads_file[1]} | \\
	(bowtie2 -q -x ${bt2_idx_path}/${params.bt2_idx_name} -U - --un-gz "${sample_id}_bt2_microbe_pangenome_unaligned.fastq.gz" \\
	-p $task.cpus --very-sensitive) 2>"${sample_id}_bt2.log" | \\
	samtools view -bS - > "${sample_id}_bt2_microbe_pangenome_aligned.bam"

	bt2_helper.sh "${sample_id}"
	
	"""
	}
}
