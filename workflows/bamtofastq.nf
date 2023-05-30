// WORKFLOW file for getting flagstats and unmapped (non-host) reads in *fastq.gz format from *.bam/*.cram files (assuming alignment to host genome)

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
    Output arguments:
      --outdir                      The output directory where the results will be saved [Default: ./pipeline_results]
      --tracedir                    The directory where nextflow logs will be saved [Default: ./pipeline_results/pipeline_info]
    AWSBatch arguments:
      --awsregion                   The AWS Region for your AWS Batch job to run on [Default: false]
      --awsqueue                    The AWS queue for your AWS Batch job to run on [Default: false]
    Others:
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


/*
========================================================================================
    Define channels for bam or cram files
========================================================================================
*/
//Just an example:
//params.rna_reads = "$baseDir/data/raw_fastq/rna/*{1,2}.{sam,bam,cram}"
//params.dna_reads = "$baseDir/data/raw_fastq/dna/*{1,2}.{sam,bam,cram}"
//https://nextflow-io.github.io/patterns/process-per-csv-record/
// The channel should look like this: [SRR493366, [/my/data/SRR493366.cram]]

if (params.process_rna && params.rna_list){
    Channel
	.fromPath( params.rna_list )
	.splitCsv(header:true) //Read in TWO column csv file with the headers: id, alignment_file
	.map { row-> tuple(row.id, tuple(file(params.rna_reads + "/" + row.alignment_file,checkIfExists: true))) }
	.set{ ch_rna_input }
} else if (params.process_rna && !params.rna_list){
	Channel.fromFilePairs( [params.rna_reads + '/**{sam,bam,cram}'], checkIfExists:true, size: 1 ).set{ ch_rna_input }
}

if (params.process_dna && params.dna_list){
    Channel
	.fromPath( params.dna_list )
	.splitCsv(header:true) //Read in TWO column csv file with the headers: id, alignment_file
	.map { row-> tuple(row.id, tuple(file(params.dna_reads + "/" + row.alignment_file,checkIfExists: true))) }
	.set{ ch_dna_input }
} else if (params.process_dna && !params.dna_list){
	Channel.fromFilePairs( [params.dna_reads + '/**{sam,bam,cram}'], checkIfExists:true, size: 1 ).set{ ch_dna_input }
}


/*
========================================================================================
    Include modules
========================================================================================
*/

include { BAM_TO_FASTQ_RNA } from '../modules/bamtofastq_rna.nf'
include { BAM_TO_FASTQ_DNA } from '../modules/bamtofastq_dna.nf'

/*
========================================================================================
    Named workflow
========================================================================================
*/

workflow BAMTOFASTQ {

if (params.process_rna){
    BAM_TO_FASTQ_RNA(rna_reads)
}
if (params.process_dna){
    BAM_TO_FASTQ_DNA(dna_reads)
}

}
