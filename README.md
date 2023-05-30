# Nextflow microbial pangenome mapping pipeline for metagenomics & metatranscriptomics 

## Introduction

Chiamh/alignment-nf is a bioinformatics pipeline that takes metagenomic **(MGX)** and/or metatranscriptomic **(MTX)** reads and aligns them to user provided pangenomes using either Bowtie2 (for bacteria) or STAR or Salmon (for eukaryotes)

The pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) , a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. 
It uses Docker containers (also compatible with Singularity) for ease of installation and computational reproducibility. 

This pipeline currently only accepts paired-end reads as inputs.

Every step in this pipeline is optional and can be toggled on or off.

## Pipeline summary for metagenomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Removal of host (human) reads by mapping to a reference genome using bwa (0.7.17) 
3. Mapping DNA reads to a microbial pangenome of choice using bowtie2 (2.4.4) followed by summaries of **unpaired** read counts per feature

## Pipeline summary for metatranscriptomic reads
1. Adapter trimming and quality control using fastp (0.22.0)
2. Removal of host (human) reads using STAR (2.7.9a), a splice aware aligner.
3. Computational removal of prokaryotic and eukaryotic rRNAs using a k-mer based strategy with bbmap (38.93)
4. Sequence de-duplication using bbmap (38.93) clumpify.sh
5. Mapping RNA reads to a microbial pangenome of choice using either bowtie2 (2.4.4) or STAR (2.7.9a) or pseudoalignment to a [decoy aware](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/) multi-species/single species transcriptome with Salmon (v 1.10.1)


## Input requirements
Either:
1. Absolute path to the **folder** containing the DNA and/or RNA reads specified with the --dna_reads and/or --rna_reads arguments. 
* This will search the folder(s) recursively for fastq files and run the pipeline on all of them.  
or:  
2. Absolute path to the **folder** containing the DNA and/or RNA reads specified with the --dna_reads and/or --rna_reads arguments **and** csv files specified with the --rna_list and --dna_list arguments.
* The csv files are in a 3 column format with headers. They correspond to the library ID, file name of read 1 and file name of read 2 respectively.
<img src='/docs/input_csv_example.PNG' width='500'>
* This will run the pipeline on the files specified in the --rna_list and/or --dna_list only 

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`) and add the nextflow executable to your $PATH

2. Install [`Docker`](https://docs.docker.com/engine/installation/)   

3. Clone the pipeline and refer to the help message
	```
	git clone https://github.com/Chiamh/alignment-nf
	
	nextflow run ./alignment-nf/main.nf --help
	```
* Add a custom config file which contains the paths to various pre-installed databases. Refer to the awsbatch.config file in this repo for an example. 
* Add a custom profile in the nextflow.config file, allowing you to specify the use of docker or singularity, and/or a task scheduler.  

4. Make sure all helper scripts in alignment-nf/bin have execute permissions

	```
	chmod +x ./alignment-nf/bin/*
	```

5. Run the pipeline

	```
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results
	```
	
* You can specify multiple profiles separated by comma, e.g. -profile docker,sge.
* Delete the work/ directory after running the pipeline to free up space taken up by intermediate files
* You can choose whether or not to save the intermediate files to the output directory with --save_intermediates (default: true)
* You have the flexibility to skip any of the steps in this pipeline. Look at the help message for more details
* If using awsbatch, you need to specify an S3 bucket as a work directory with the -bucket-dir argument

Further usage

	```
	
	#Runs the pipeline but skips preprocessing with fastp, host removal (decontamination), rRNA removal and RNA sequence de-duplication.
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results \
	--preprocess false --decont false --remove_rRNA false --dedupe false --map true --rna_mapper bowtie2 --save_intermediates false
	
	#Runs the pipeline but only process DNA reads (metagenomes)
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --dna_reads /path/to/metagenomes --outdir /path/to/results --process_rna false
	
	#Runs the pipeline but only process RNA reads (metatranscriptomes)
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --outdir /path/to/results --process_dna false
	
	#Runs the pipeline using STAR to map RNA reads to a eukaryotic/fungal pangenome.
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --outdir /path/to/results --process_dna false --rna_mapper star
	
	#Runs the pipeline using Salmon to map RNA reads to a eukaryotic transcriptome (can be multi-species).
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile --rna_reads /path/to/metatranscriptomes --outdir /path/to/results --process_dna false --rna_mapper salmon
	
	```

## Output files

* **Caveat: RNA seq mapping with bowtie2 uses inputs in a "non-paired" fashion despite read pairs**
* **Caveat: RNA seq mapping with STAR uses inputs in a "paired" fashion**
* **Caveat: DNA seq mapping ALWAYS uses inputs in a "non-paired" fashion despite read pairs**

It is preferable to map metagenomic reads to pangenomes (bacterial or eukaryotic) in a **non-paired** fashion despite paired-end data due to variable gene order in microbial strains.

It is preferable to map reads in a **non-paired** fashion to bacterial pangenomes despite paired-end data due to polycistronic RNAs. [Read this.](https://github.com/biobakery/humann#humann-30-and-paired-end-sequencing-data)

Mapping **RNA reads with STAR** assumes a eukaryotic/fungal genome or pangenome, mono-cistronic RNAs and will **use paired reads and output paired read counts**. 

Pseudoalignment of **RNA reads with Salmon** to multi-species transcriptome is also possible. This is especially useful for eukaryotic species without a well-annotated pangenome.

* reads/DNA (for metagenomes) or reads/RNA (for metatranscriptomes)
    * These folders contain processed reads in fastq.gz format. "Processed" means adapter removal, host read removal, rRNA removal (for MTX) and/or de-duplication (for MTX)
	* If --save_intermediates is true, there will be a temp/ subfolder containing the intermediate files
	* \*_fastp_{1,2}.fastq.gz are reads that have been preprocessed by fastp. This step occurs before decontamination.
	* \*_host_unmapped_{1,2}.fastq.gz are reads that do not map to host reference genome. The removal of host reads is called "decontamination" in this workflow.
	* \*_mRNA_{1,2}.fastq.gz are reads where human, fungal and bacterial rRNAs have been computationally removed. In the full workflow, rRNA removal occurs after decontamination.
	* \*_dedup_{1,2}.fastq.gz are reads that have been de-duplicated. In the full workflow, this step occurs after rRNA removal. Only applicable for RNA seq.

* bt2_out/DNA (for metagenomes) or bt2_out/RNA (for metatranscriptomes) for pangenome mapping results using bowtie2
    * \*bt2_microbe_pangenome_aligned.bam : BAM file after alignment of reads (non-paired manner) to pangene catalog.
	* \*bt2_microbe_pangenome_filtered_cov.tsv : Tab separated file containing unpaired read coverage across pangenes. Only pangenes with >= 50% coverage are reported here.
	* \*bt2_microbe_pangenome.fastq.gz : All reads that did not align to the pangene catalog.

* STAR_out/ (for metatranscriptomes only) RNA aligned to pangenome using STAR
    * \*star_microbe_pangenome_aligned.bam : BAM file after alignment of reads (PAIRED manner) to pangene catalog.
	* \*star_microbe_pangenome_filtered_cov.tsv : Tab separated file containing PAIRED read coverage across pangenes. Only pangenes with >= 50% coverage are reported here.
	* \*star_microbe_pangenome.fastq.gz : All reads that did not align to the pangene catalog.
	
* salmon_out/ (for metatranscriptomes only) RNA pseudoalignment to transcriptome using Salmon
	* \*_quant.sf 

## Alternative workflows

### Concatenate
The "concatenate" subworkflow merges fastq.gz files across different lanes for the same sample ID.\
**Inputs**\
Absolute path to the **folder** containing the DNA and/or RNA reads specified with the --dna_reads and/or --rna_reads arguments. 
* This will search the folder(s) recursively for fastq files and run the pipeline on all of them.

	```
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile -entry concatenate --rna_reads /path/to/metatranscriptomes --dna_reads /path/to/metagenomes --outdir /path/to/results
	
	```
	
### bamtofastq
The "bamtofastq" subworkflow extracts non-host reads (both R1 and R2 unmapped) from sam/bam/cram files which were produced by alignment to a host reference.
* Outputs fastq.gz files.
* Also generates flagstat.tsv files. The first column contains the values for QC-passed reads (not read pairs!), the second column has the values for QC-failed reads and the third contains the category names.
* Singleton non-host reads are saved to a separate results folder.
* Important: This subworkflow needs the user to specify 'input_is_bam true' in the command\
	
**Inputs**\
Either:  
1. Absolute path to the **folder** containing the DNA and/or RNA reads specified with the --dna_reads and/or --rna_reads arguments. 
* This will search the folder(s) recursively for sam/bam/cram files and run the pipeline on all of them.  
or:  
2. Absolute path to the **folder** containing the DNA and/or RNA reads specified with the --dna_reads and/or --rna_reads arguments **and** csv files specified with the --rna_list and --dna_list arguments.
* The csv files are in a two column format with headers. They correspond to the library ID, and file name of alignment file.
* This will run the pipeline on the files specified in the --rna_list and/or --dna_list only 
<img src='/docs/input_csv_example2.PNG' width='183' height='63'>

	```
	nextflow run ./alignment-nf/main.nf -profile docker,your_profile -entry bamtofastq --input_is_bam true --rna_reads /path/to/RNA_bam_files --dna_reads /path/to/DNA_bam_files --outdir /path/to/results
	
	```
## Contact

Minghao Chia: chia_minghao@gis.a-star.edu.sg, chiaminghao@gmail.com
