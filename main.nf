#!/usr/bin/env nextflow

/*
========================================================================================
    Meta-omic sequence alignment
========================================================================================
    Github : https://github.com/Chiamh/alignment-nf
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

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
    
    The main workflow can take up a lot of disk space with intermediate fastq files. 
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


/*
========================================================================================
    Include modules
========================================================================================
*/

include { FULL } from './workflows/full_workflow.nf'

/*
========================================================================================
    Main workflow (default)
========================================================================================
*/

// this is the main workflow
workflow {
    
    FULL ()
     
}