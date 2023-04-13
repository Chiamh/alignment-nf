// Reads are mapped as paired end because eukaryotic reads are not poly-cistronic.

process SALMON {
	label "process_medium"
	tag "${sample_id}"
	publishDir "${params.outdir}/salmon_out", mode: 'copy'
	
	input:
	path salmon_index
	tuple val(sample_id), path(reads_file)
	
	output:
	tuple val(sample_id), path("${sample_id}"), emit: results_folder
	tuple val(sample_id), path("*_quant.sf"), emit: results
	tuple val(sample_id), path("*info.json"), emit: json_info
	
	when:
	params.map && params.process_rna
	
	script:

	"""
		
	salmon quant -i ${salmon_index} -l A -1 ${reads_file[0]} -2 ${reads_file[1]} --validateMappings -o ${sample_id}
	
	if [ -f ${sample_id}/quant.sf ]; then
        cp ${sample_id}/quant.sf "${sample_id}_quant.sf"
    fi
	
	
    if [ -f ${sample_id}/aux_info/meta_info.json ]; then
        cp ${sample_id}/aux_info/meta_info.json "${sample_id}_meta_info.json"
    fi
	
			
	"""
}

