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
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main' 
include { SEQTK_TRIM             } from '../modules/nf-core/seqtk/trim/main' 
include { STAR_ALIGN } from '../modules/nf-core/star/align/main'
include { HISAT2_BUILD } from '../modules/nf-core/hisat2/build/main' 
include { HISAT2_ALIGN } from '../modules/nf-core/hisat2/align/main'  
include { HISAT2_EXTRACTSPLICESITES } from '../modules/nf-core/hisat2/extractsplicesites/main' 
include { STAR_GENOMEGENERATE } from '../modules/nf-core/star/genomegenerate/main'                                                          
include { GUNZIP } from '../modules/nf-core/gunzip/main'
<<<<<<< Updated upstream
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main' 
=======
include { getGenomeAttribute } from '../subworkflows/local/utils_nfcore_rnasequencing_pipeline'

>>>>>>> Stashed changes
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQUENCING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
<<<<<<< Updated upstream
    ch_fasta       // string: path to fasta file from params
    ch_gtf         // string: path to gtf file from params
=======
>>>>>>> Stashed changes

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    FASTQC (
        ch_samplesheet
    )

<<<<<<< Updated upstream
=======
    files_to_align = ch_samplesheet

>>>>>>> Stashed changes
    if (params.trimmer == "none") {
        println "No trimming will be performed"
    } else {
        if (params.trimmer == "seqtk_trim") {
            println "Using seqtk for trimming"
            SEQTK_TRIM (
                ch_samplesheet
            )
<<<<<<< Updated upstream
            FASTQC_AFTER_TRIM(
                SEQTK_TRIM.out.reads
            )
            ch_versions.mix(SEQTK_TRIM.out.versions)
=======
            //ch_multiqc_files = ch_multiqc_files.mix(trimmed_files)
            ch_versions = ch_versions.mix(seqtk_versions)
>>>>>>> Stashed changes
        } else if (params.trimmer == "trimgalore") {
            println "Using TrimGalore for trimming"
            TRIMGALORE (
                ch_samplesheet
            )
<<<<<<< Updated upstream
            FASTQC_AFTER_TRIM(
                TRIMGALORE.out.reads
            )
        }
    }

    ch_fasta_for_index = Channel.value([['id':'genome'], file(ch_fasta, checkIfExists: true)])
    ch_gtf_for_index = Channel.value([['id':'gtf'], file(ch_gtf, checkIfExists: true)])

=======

        }

        
        ch_after_trim_fastqc_reports = FASTQC_AFTER_TRIM (
            trimmed_files 
        )
    }
    ch_fasta_for_index = Channel.value([['id':'genome'], file(getGenomeAttribute('fasta'), checkIfExists: true)])
    ch_gtf_for_index = Channel.value([['id':'gtf'], file(getGenomeAttribute('gtf'), checkIfExists: true)])
>>>>>>> Stashed changes
    if (params.aligner == 'star') {
        // STAR
        println "Using STAR for alignment"

        STAR_GENOMEGENERATE (
            ch_fasta_for_index,
            ch_gtf_for_index
        )
<<<<<<< Updated upstream

        GUNZIP (
            ch_samplesheet.map { it[1] }
=======
        ch_samplesheet.view()
        ch_paths = ch_samplesheet.map { it[1] }
    
        (a,b) = GUNZIP (
            ch_paths
>>>>>>> Stashed changes
        )

        STAR_ALIGN (
            GUNZIP.out.gunzip,
            STAR_GENOMEGENERATE.out.index,
            ch_gtf_for_index,
            false,
            "Illumina",
            "Dummy"
<<<<<<< Updated upstream
=======
        )}
    else if (params.aligner == 'hisat2') {
        //HISAT2
        ch_gtf = Channel.value([['id':'genome'], file(ch_gtf, checkIfExists: true)])
        println "Type of ch_gtf: ${ch_gtf.getClass()}"
        println "Type of params.gtf: ${params.gtf.getClass()}"

        (txt, versions) = HISAT2_EXTRACTSPLICESITES(
            ch_gtf
>>>>>>> Stashed changes
        )

        SAMTOOLS_FAIDX (
            ch_fasta_for_index,
            [[],[]], //Channel.of([]),
            false
        )

        SAMTOOLS_SORT(
            STAR_ALIGN.out.bam,
            ch_fasta_for_index,
            "bai"
        )

        PICARD_MARKDUPLICATES (
            SAMTOOLS_SORT.out.bam,
            ch_fasta_for_index,
            SAMTOOLS_FAIDX.out.fai
        )
    } else if (params.aligner == 'hisat2') {
        //HISAT2
        HISAT2_EXTRACTSPLICESITES(
            ch_gtf_for_index
        )

        HISAT2_BUILD (
            ch_fasta_for_index,
            ch_gtf_for_index,
            HISAT2_EXTRACTSPLICESITES.out.txt
        )

        HISAT2_ALIGN (
            ch_samplesheet,
            HISAT2_BUILD.out.index,
            HISAT2_EXTRACTSPLICESITES.out.txt
        )
    } else {
        println "Please select a valid aligner: 'star' or 'hisat2'"
    }


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
