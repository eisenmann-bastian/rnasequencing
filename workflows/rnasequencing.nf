/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER_TRIM } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'      
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnasequencing_pipeline'

include { SEQTK_TRIM             } from '../modules/nf-core/seqtk/trim/main' 
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { HISAT2_BUILD } from '../modules/nf-core/hisat2/build/main' 
include { HISAT2_ALIGN } from '../modules/nf-core/hisat2/align/main'  
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main' 
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'                                                          
include { GUNZIP } from '../modules/nf-core/gunzip/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQUENCING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_trimmer     // string: trimmer choice from params
    ch_fasta       // string: path to fasta file from params
    ch_gtf         // string: path to gtf file from params

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_pre_trim_fastqc_reports = FASTQC (
        ch_samplesheet
    )

    files_to_align = ch_samplesheet

    if (ch_trimmer == "none") {
        println "No trimming will be performed"
        files_to_align = ch_samplesheet
    } else {
        if (ch_trimmer == "seqtk_trim") {
            println "Using seqtk for trimming"
            (trimmed_files, seqtk_versions) = SEQTK_TRIM (
                ch_samplesheet
            )
            //ch_multiqc_files = ch_multiqc_files.mix(trimmed_files)
            ch_versions = ch_versions.mix(seqtk_versions)
        } else if (ch_trimmer == "trimgalore") {
            println "Using TrimGalore for trimming"
            (trimmed_files, logs, unpaired, html, zip, versions) = TRIMGALORE (
                ch_samplesheet
            )
            /*trimmed_files.view()
            ch_multiqc_files = ch_multiqc_files.mix(trimmed_files)
            ch_versions = ch_versions.mix(versions)*/
        }

        
        ch_after_trim_fastqc_reports = FASTQC_AFTER_TRIM (
            trimmed_files
        )
    }
    
    ch_fasta_for_index = Channel.value([['id':'genome'], file(ch_fasta, checkIfExists: true)])
    ch_gtf_for_index = Channel.value([['id':'gtf'], file(ch_gtf, checkIfExists: true)])

    // memory issues  
    (index,versions) = STAR_GENOMEGENERATE (
        ch_fasta_for_index,
        ch_gtf_for_index
    )
    /*ch_star_index = STAR_GENOMEGENERATE.out.index
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

    (index,versions)=STAR_GENOMEGENERATE {
        ch_star_refs,
        ch_star_gtfs
    }*/
    ch_samplesheet.view()
    ch_paths = ch_samplesheet.map { it[1] }
   
    (a,b) = GUNZIP (
        ch_paths
    )
    a.view()
    b.view()
    STAR_ALIGN (
        a,
        index,
        ch_gtf_for_index,
        false,
        "Illumina",
        "Dummy"
    )


    /* HISAT2
    ch_gtf = Channel.value([['id':'genome'], file(ch_gtf, checkIfExists: true)])

    (txt, versions) = HISAT2_EXTRACTSPLICESITES(
        ch_gtf
    )

    (index, versions) = HISAT2_BUILD (
        ch_fasta_for_index,
        ch_gtf,
        txt
    )

    HISAT2_ALIGN (
        files_to_align,
        index,
        txt
    )*/
   
    

    /*
    (trimmed_files,seqtk_versions) = SEQTK_TRIM (
        ch_samplesheet
    )
    trimmed_files.view()

    
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    */

   

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnasequencing_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
   
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
