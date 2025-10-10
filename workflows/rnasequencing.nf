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
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_SORT } from '../modules/nf-core/samtools/sort/main' 

include { getGenomeAttribute } from '../subworkflows/local/utils_nfcore_rnasequencing_pipeline'
include { SUBREAD_FEATURECOUNTS } from '../modules/nf-core/subread/featurecounts/main'
include { GFFREAD } from '../modules/nf-core/gffread/main'
include { SALMON_INDEX } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT } from '../modules/nf-core/salmon/quant/main'                                                

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQUENCING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_fasta_for_index = Channel.value([['id':'genome'], file(getGenomeAttribute('fasta'), checkIfExists: true)])
    ch_gtf_for_index = Channel.value([['id':'gtf'], file(getGenomeAttribute('gtf'), checkIfExists: true)])

    if (params.run_fastqc_at_start) {
        println "Running FastQC on raw reads"
        
        FASTQC (
            ch_samplesheet
        )
    }

    if (params.trimmer == "none") {
        println "No trimming will be performed"

        ch_reads_to_align = ch_samplesheet
    } else {
        if (params.trimmer == "seqtk_trim") {
            println "Using seqtk for trimming"

            SEQTK_TRIM (
                ch_samplesheet
            )
            
            ch_trimmed_reads = SEQTK_TRIM.out.reads
        } else if (params.trimmer == "trimgalore") {
            println "Using TrimGalore for trimming"

            TRIMGALORE (
                ch_samplesheet
            )

            ch_trimmed_reads = TRIMGALORE.out.reads
        }
        if (params.run_fastqc_after_trim) {
            println "Running FastQC on trimmed reads"

            FASTQC_AFTER_TRIM(
                ch_trimmed_reads
            )
        }

        ch_reads_to_align = ch_trimmed_reads
    }

    //
    // PSEUDO-ALIGNER: Transcript quantification with SALMON
    //
    if (params.pseudo_aligner == 'salmon') {
        println "Running SALMON for transcript quantification"
        
        // Generate transcriptome FASTA from genome and GTF using GFFREAD
        GFFREAD (
            ch_gtf_for_index.map{meta, gtf -> [meta, gtf]},    // tuple val(meta), path(gff)
            ch_fasta_for_index.map{it[1]}                      // path fasta
        )
        ch_versions = ch_versions.mix(GFFREAD.out.versions.first())
        
        // Create SALMON index with transcriptome as main target, genome as decoy
        SALMON_INDEX (
            ch_fasta_for_index.map{it[1]},                     // genome_fasta (as decoy)
            GFFREAD.out.gffread_fasta.map{it[1]}               // transcript_fasta (from GFFREAD)
        )
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())

        // Quantify transcripts with SALMON
        SALMON_QUANT (
            ch_reads_to_align,                                 // tuple val(meta), path(reads)
            SALMON_INDEX.out.index,                            // path index
            ch_gtf_for_index.map{it[1]},                       // path gtf
            GFFREAD.out.gffread_fasta.map{it[1]},              // path transcript_fasta (from GFFREAD)
            false,                                             // val alignment_mode
            []                                                 // val lib_type (auto-detect)
        )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]})
    }

    if (params.aligner == 'star') {
        // STAR
        println "Using STAR for alignment"

        STAR_GENOMEGENERATE (
            ch_fasta_for_index,
            ch_gtf_for_index
        )

        GUNZIP (
            ch_reads_to_align.map { it[1] }
        )

        STAR_ALIGN (
            GUNZIP.out.gunzip,
            STAR_GENOMEGENERATE.out.index,
            ch_gtf_for_index,
            false,
            "Illumina",
            "Dummy"
        )

        SAMTOOLS_FAIDX (
            ch_fasta_for_index,
            [[],[]],
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
        HISAT2_EXTRACTSPLICESITES(
            ch_gtf_for_index
        )

        HISAT2_BUILD (
            ch_fasta_for_index,
            ch_gtf_for_index,
            HISAT2_EXTRACTSPLICESITES.out.txt
        )

        HISAT2_ALIGN (
            ch_reads_to_align,
            HISAT2_BUILD.out.index,
            HISAT2_EXTRACTSPLICESITES.out.txt
        )

        ch_bam = HISAT2_ALIGN.out.bam.map { it[1] }
        ch_gtf = ch_gtf_for_index.map { it[1] }

        ch_featurecounts_input = ch_bam.combine(ch_gtf).map { bam, gtf ->
            tuple([ id: 'genome' ], bam, gtf)
        }

        SUBREAD_FEATURECOUNTS(ch_featurecounts_input)                                             
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
