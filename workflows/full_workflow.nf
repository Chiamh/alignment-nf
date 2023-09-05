/*
========================================================================================
    Help messages and warnings
========================================================================================
*/

//https://www.baeldung.com/groovy-def-keyword
//https://www.nextflow.io/blog/2020/cli-docs-release.html

def helpMessage() {
  // adapted from nf-core
    log.info"""
    
    Usage for main workflow:
    The typical command for running the pipeline is as follows:
      nextflow run main.nf
	  --rna_list PATH_TO_INPUT_CSV_FOR_RNA
	  --dna_list PATH_TO_INPUT_CSV_FOR_DNA
	  --rna_reads FOLDER_FOR_RNA_READS
      --dna_reads FOLDER_FOR_DNA_READS
	  --rna_mapper CHOICE_OF_MAPPER_FOR_RNASEQ_ALIGNMENT
    
	NOTE: A more user-friendly approach is to specify these parameters in a *.config file under a custom profile 
	
    IMPT: Set either the --process_rna or --process_dna arguments to false if no RNA or DNA reads are provided, respsectively. 
    
    The main workflow can take up a lot of disk space with intermediate fastq files. Be sure to clean up the work directory after the pipeline has finished.
	You can set the --save_intermediates flag to false to avoid copying large intermediates in the output folder. 
           
    Input and database arguments are null by default.
    Rather than manually specifying the paths to so many databases, it is best to create a custom nextflow config file.
     
    Input arguments:
      --rna_list                    Path to a three column csv file with headers: id,read1,read2 for metatranscriptomic reads. If not defined, workflow will search input folder for all valid input fastq files. 
	  --dna_list                    Path to a three column csv file with headers: id,read1,read2 for metagenomic reads. If not defined, workflow will search input folder for all valid input fastq files.
	  --rna_reads                   Path to a folder containing all input metatranscriptomic reads (this will be recursively searched for *fastq.gz/*fq.gz/*fq/*fastq files)
      --dna_reads                   Path to a folder containing all input metagenomic reads (this will be recursively searched for *fastq.gz/*fq.gz/*fq/*fastq files)
    Database arguments:
      --bwaidx_path                 Path to the folder with host (e.g. human) reference genome and bwa index for decont
      --bwaidx			            Name of the bwa index for decont e.g. hg38.fa
	  --decont_star_index           Path to the directory containing the index for the host (e.g. human) genome for STAR aligner
      --star_index                  Path to the directory containing the index for pangenome mapping using STAR
      --ribokmers                   Path to the eukaryotic and prokaryotic ribokmer database for computational rRNA removal using BBmap
      --bt2_idx_path                Path to the folder with bowtie2 index for custom-built microbial pangenome/gene catalog
      --bt2_idx_name                Name of the bowtie2 index for the pangenome/gene catalog e.g. IHSMGC
	  --salmon_index                Path to the folder with the salmon index for decoy aware multi-species or single species transcriptome
    Workflow options:
      --process_rna                 Turns on steps to process metatranscriptomes [Default: true]. If true, --rna_reads is a mandatory argument
      --process_dna                 Turns on steps to process metagenomes [Default: true]. If true, --dna_reads is a mandatory argument
      --preprocess                  Preprocess reads with fastp (trimming and QC) [Default: true]
	  --decont						Map reads to reference genomes e.g. host/human to remove them [Default: true]
	  --remove_rRNA                 Remove ribosomal RNA reads
	  --dedupe                      Deduplication for metatranscriptomes [Default: true]
	  --map                         Map MGX and MTX reads to pangenomes [Default: true]
      --rna_mapper                  Choice of which mapper to use for metatranscriptome read mapping. Choose from: bowtie2, star or salmon [Default: bowtie2]
	  --save_intermediates          Copy intermediate files to output directory [Default: true]
      --force_paired                forces bowtie2 to map reads in paired end mode for both mtx and mgx [Default: false]
    Output arguments:
      --outdir                      The output directory where the results will be saved [Default: ./pipeline_results]
      --tracedir                    The directory where nextflow logs will be saved [Default: ./pipeline_results/pipeline_info]
    AWSBatch arguments:
      --awsregion                   The AWS Region for your AWS Batch job to run on [Default: false]
      --awsqueue                    The AWS queue for your AWS Batch job to run on [Default: false]
    Others:
      --cleanup			    This option will enable nextflow work/temp folder cleanup upon pipeline completion (no errors).
                                    All intermediate files from nexftlow processes' workdirs will be cleared. 
				    It will not clear cached folders coming from previous pipeline runs. [Default: false]
      --help		            Display this help message
    """
}

if (params.help){
    helpMessage()
    exit 0
}



if (!params.rna_reads && !params.dna_reads){
    helpMessage()
    log.info"""
    [Error] The path to at least one input folder for sequences is required
    """.stripIndent()
    exit 0
}

if (!params.rna_reads && params.process_rna){
    helpMessage()
    log.info"""
    [Error] The path to input RNA sequences is required because --process_rna is true
    """.stripIndent()
    exit 0
}

if (!params.dna_reads && params.process_dna){
    helpMessage()
    log.info"""
    [Error] The path to input DNA sequences is required because --process_dna is true
    """.stripIndent()
    exit 0
}


if (!params.bwaidx_path && params.decont && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --bwaidx_path is required for removal of host (human) reads from metagenomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.bwaidx && params.decont && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --bwaidx is required for removal of host (human) reads from metagenomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}


if (!params.decont_star_index && params.decont && params.process_rna){
    helpMessage()
    log.info"""
    [Error] --decont_star_index is required for removal of host (e.g. human) reads from metatranscriptomic sequences (decontamination step)
    """.stripIndent()
    exit 0
}

if (!params.ribokmers && params.remove_rRNA && params.process_rna){
    helpMessage()
    log.info"""
    [Error] --ribokmers is required for removal of rRNA reads
    """.stripIndent()
    exit 0
}

if (params.map && params.process_rna && params.rna_mapper == null ){
    helpMessage()
    log.info"""
    [Error] Specify a valid --rna_mapper either bowtie2, star or salmon
    """.stripIndent()
    exit 0
}

if (params.map && params.process_rna && params.rna_mapper != 'bowtie2' && params.rna_mapper != 'star' && params.rna_mapper != 'salmon'){
    helpMessage()
    log.info"""
    [Error] Specify a valid --rna_mapper either bowtie2, star or salmon
    """.stripIndent()
    exit 0
}


if (params.map && !params.bt2_idx_path && params.process_rna && params.rna_mapper == 'bowtie2'){
    helpMessage()
    log.info"""
    [Error] --bt2_idx_path is required for mapping of MTX reads to pangenome with bowtie2
    """.stripIndent()
    exit 0
}

if (params.map && !params.bt2_idx_name && params.process_rna && params.rna_mapper == 'bowtie2'){
    helpMessage()
    log.info"""
    [Error] --bt2_idx_name is required for mapping of MTX reads to pangenome with bowtie2
    """.stripIndent()
    exit 0
}

if (params.map && !params.star_index && params.process_rna && params.rna_mapper == 'star'){
    helpMessage()
    log.info"""
    [Error] --star_index is required for mapping of MTX reads to gene catalog with STAR
    """.stripIndent()
    exit 0
}

if (params.map && !params.salmon_index && params.process_rna && params.rna_mapper == 'salmon'){
    helpMessage()
    log.info"""
    [Error] --salmon_index is required for pseudoalignment of MTX reads to transcriptome
    """.stripIndent()
    exit 0
}

if (params.map && !params.bt2_idx_path && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --bt2_idx_path is required for mapping of MGX reads to pangenome with bowtie2
    """.stripIndent()
    exit 0
}

if (params.map && !params.bt2_idx_name && params.process_dna){
    helpMessage()
    log.info"""
    [Error] --bt2_idx_name is required for mapping of MGX reads to pangenome with bowtie2
    """.stripIndent()
    exit 0
}



/*
========================================================================================
    Define channels for read pairs
========================================================================================
*/
//Just an example:
//params.rna_reads = "$baseDir/data/raw_fastq/rna/*{1,2}.{fq.gz,fastq.gz}"
//params.dna_reads = "$baseDir/data/raw_fastq/dna/*{1,2}.{fq.gz,fastq.gz}"
//https://nextflow-io.github.io/patterns/process-per-csv-record/
// The channel should look like this: [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]

if (params.process_rna && params.rna_list && !params.input_is_bam){
    Channel
	.fromPath( params.rna_list )
	.splitCsv(header:true) //Read in 3 column csv file with the headers: id, read1 and read2
	.map { row-> tuple(row.id, tuple(file(params.rna_reads + "/" + row.read1,checkIfExists: true), file(params.rna_reads + "/" + row.read2,checkIfExists: true))) }
	.set{ ch_rna_input }
} else if (params.process_rna && !params.rna_list && !params.input_is_bam){
	Channel.fromFilePairs( [params.rna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_rna_input }
}

if (params.process_dna && params.dna_list && !params.input_is_bam){
    Channel
	.fromPath( params.dna_list )
	.splitCsv(header:true) //Read in 3 column csv file with the headers: id, read1 and read2
	.map { row-> tuple(row.id, tuple(file(params.dna_reads + "/" + row.read1,checkIfExists: true), file(params.dna_reads + "/" + row.read2,checkIfExists: true))) }
	.set{ ch_dna_input }
} else if (params.process_dna && !params.dna_list && !params.input_is_bam){
	Channel.fromFilePairs( [params.dna_reads + '/**{R,.,_}{1,2}*{fastq,fastq.gz,fq,fq.gz}'], checkIfExists:true ).set{ ch_dna_input }
}


/*
========================================================================================
    Include modules
========================================================================================
*/

include { FASTP_DECONT_RNA } from '../modules/fastp_decont_rna.nf'
include { FASTP_DECONT_DNA } from '../modules/fastp_decont_dna.nf'
include { FASTP_RNA } from '../modules/fastp_rna.nf'
include { FASTP_DNA } from '../modules/fastp_dna.nf'
include { RIBOFILTER } from '../modules/rRNAfilter.nf'
include { DEDUP } from '../modules/dedup.nf'
include { BT2ALIGN_RNA } from '../modules/bt2align_rna.nf'
include { BT2ALIGN_DNA } from '../modules/bt2align_dna.nf'
include { STARALIGN_RNA } from '../modules/staralign_rna.nf'
include { SALMON } from '../modules/salmon.nf'
include { DECONT_DNA } from '../modules/decont_dna.nf'
include { DECONT_RNA } from '../modules/decont_rna.nf'

/*
========================================================================================
    Workflow
========================================================================================
*/

//https://bioinformatics.stackexchange.com/questions/15739/use-conditional-in-workflow-in-nextflow-dsl2
//https://github.com/nextflow-io/patterns/blob/master/docs/conditional-process.adoc

// this is the main workflow
workflow FULL {
    if ( params.process_rna ){
	
		if ( params.preprocess && params.decont ){
			FASTP_DECONT_RNA(params.decont_star_index, ch_rna_input)
			ch_qc_reads = FASTP_DECONT_RNA.out.reads
			} else if ( params.preprocess && !params.decont ){
				FASTP_RNA(ch_rna_input)
				ch_qc_reads = FASTP_RNA.out.reads
				} else if ( !params.preprocess && params.decont ){
					DECONT_RNA(params.decont_star_index, ch_rna_input)
					ch_qc_reads = DECONT_RNA.out.reads
					} else if ( !params.preprocess && !params.decont ){
						ch_qc_reads = ch_rna_input
						}
    
		if ( params.remove_rRNA ){
			RIBOFILTER(params.ribokmers, ch_qc_reads)
		} 
	
		if ( params.remove_rRNA && params.dedupe ){
			DEDUP(RIBOFILTER.out.reads)
			ch_rna_decont = DEDUP.out.reads
			} else if ( !params.remove_rRNA && params.dedupe ){
				DEDUP(ch_qc_reads)
				ch_rna_decont = DEDUP.out.reads
				} else if ( !params.remove_rRNA && !params.dedupe ){
					ch_rna_decont = ch_qc_reads
					} else if ( params.remove_rRNA && !params.dedupe ){
						ch_rna_decont = RIBOFILTER.out.reads
						}
	
		if ( params.map && params.rna_mapper == 'bowtie2' ){
			BT2ALIGN_RNA(params.bt2_idx_path, ch_rna_decont)
			} else if ( params.map && params.rna_mapper == 'star' ){
				STARALIGN_RNA(params.star_index, ch_rna_decont)
				} else if ( params.map && params.rna_mapper == 'salmon' ) { 
					SALMON(params.salmon_index, ch_rna_decont)
				}
	}
    
    if ( params.process_dna ){
        
		if ( params.preprocess && params.decont){
			FASTP_DECONT_DNA(params.bwaidx_path, ch_dna_input)
			ch_dna_decont = FASTP_DECONT_DNA.out.reads
		} else if ( params.preprocess && !params.decont ){
				FASTP_DNA(ch_dna_input)
				ch_dna_decont = FASTP_DNA.out.reads
			} else if ( !params.preprocess && params.decont ){
					DECONT_DNA(params.bwaidx_path, ch_dna_input)
					ch_dna_decont = DECONT_DNA.out.reads
				} else if ( !params.preprocess && !params.decont ){
					ch_dna_decont = ch_dna_input
					} 
		
		if ( params.map ){
			BT2ALIGN_DNA(params.bt2_idx_path, ch_dna_decont)
		}
	}
}
