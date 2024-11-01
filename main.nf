#!/usr/bin/env Nextflow

nextflow.enable.dsl = 2

def prepare_input (LinkedHashMap row) {
    
    // define inputs
    accession = row.BioSample_accession
    id        = row.sample_id
    read1     = row.read1_dehosted
    read2     = row.read2_dehosted


    // check for required inputs
    if( ! accession ){ exit 1, "ERROR: 'BioSample_accession' missing for ${row}" }
    if( ! id        ){ exit 1, "ERROR: 'sample_id' missing for ${row}"           }
    if( ! read1     ){ exit 1, "ERROR: 'read1_dehosted' missing for ${row}"      }
    if( ! read2     ){ exit 1, "ERROR: 'read2_dehosted' missing for ${row}"      }

    // check that read files exist
    if( ! file(read1).exists() ){ exit 1, "ERROR: ${read1} does not exist!" }
    if( ! file(read2).exists() ){ exit 1, "ERROR: ${read2} does not exist!" }

    // check that accession and id look right
    if(  ! accession ==~ /^SAMN\d{8}$/ ){ exit 1, "ERROR: ${accession} does not look like a NCBI BioSample accession." }
    if( ! id ==~ /WA-PHL-\d{6}$/ ){ exit 1, "ERROR: ${id} does not look like a WA PHL ID." }

    // create output
    result = [ accession: accession, id: id, read1: file(read1), read2: file(read2) ]

    return result
}

process STAGE_READS {
    publishDir "${params.stagedir}"
    container 'public.ecr.aws/o8h2f0o1/bigbacter-base:1.0.0'

    input:
    tuple val(id), path(read1, stageAs: "reads/*"), path(read2, stageAs: "reads/*")

    output:
    path "*.fastq.gz"
    tuple val(id), env(SEQUENCER), emit: sequencer

    script:
    """
    # rename files
    ## read 1
    mv ${read1} ${id}_R1.fastq.gz
    ## read 2
    mv ${read2} ${id}_R2.fastq.gz

    # determine the sequencer
    case \$(zcat ${read1} | head -n1 | sed 's/:.*//g' | cut -c 2) in
        "V")
            SEQUENCER="NextSeq 2000"
            ;;
        "M")
            SEQUENCER="Illumina MiSeq"
            ;;
        *)
            SEQUENCER="Unknown"
            ;;
    esac
    """
}

process PUBLISH_META {
    publishDir "${params.outdir}", mode: 'copy'
    container 'public.ecr.aws/o8h2f0o1/bigbacter-base:1.0.0'

    input:
    path meta

    output:
    path meta

    script:
    """
    """
}

workflow {

    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ prepare_input(it) }
        .set{ ch_input }

    STAGE_READS (
        ch_input.map{ it -> [ it.id, it.read1, it.read2 ] }
    )

    ch_input
        .map{ it -> [ it.id, it ] }
        .join( STAGE_READS.out.sequencer, by: 0 )
        .map{ id, it, sequencer -> it + [sequencer: sequencer] }
        .map{ it -> "${it.accession},${it.id},Baseline surveillance (random sampling) of severe acute respiratory syndrome coronavirus 2,WGS,VIRAL RNA,PCR,paired,ILLUMINA,${it.sequencer},Whole genome sequencing (tiled-amplicon) of severe acute respiratory syndrome coronavirus 2,fastq,${it.id}_R1.fastq.gz,${it.id}_R2.fastq.gz" }
        .set{ ch_meta }
    Channel.of("biosample_accession,library_ID,title,library_strategy,library_source,library_selection,library_layout,platform,instrument_model,design_description,file_type,filename,filename2")
        .concat( ch_meta )
        .collectFile(name: "sra_metadata.csv", sort: 'index', newLine: true)
        .set{ ch_meta }
    
    PUBLISH_META (
        ch_meta
    )

}

