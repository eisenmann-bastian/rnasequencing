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

    ch_fasta = file(getGenomeAttribute('fasta'), checkIfExists: true)
    ch_gtf = file(getGenomeAttribute('gtf'), checkIfExists: true)

    ch_fasta_for_index = Channel.value([['id':'genome'], ch_fasta])
    ch_gtf_for_index = Channel.value([['id':'gtf'], ch_gtf])

    if (params.run_fastqc_at_start) {
        println "Running FastQC on raw reads"
        
        FASTQC (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
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
            ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})

            ch_trimmed_reads = TRIMGALORE.out.reads
        }
        if (params.run_fastqc_after_trim) {
            println "Running FastQC on trimmed reads"

            FASTQC_AFTER_TRIM(
                ch_trimmed_reads
            )
            ch_versions = ch_versions.mix(FASTQC_AFTER_TRIM.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_AFTER_TRIM.out.zip.collect{it[1]})
        }

        ch_reads_to_align = ch_trimmed_reads
    }

    if (params.mode == 'salmon') {
        println "Running SALMON for transcript quantification"

        GFFREAD (
            ch_gtf_for_index.map{meta, gtf -> [meta, gtf]},    
            ch_fasta_for_index.map{it[1]}                      
        )
        ch_versions = ch_versions.mix(GFFREAD.out.versions.first())
        
        SALMON_INDEX (
            ch_fasta_for_index.map{it[1]},                     
            GFFREAD.out.gffread_fasta.map{it[1]}               
        )
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions.first())

        SALMON_QUANT (
            ch_reads_to_align,                               
            SALMON_INDEX.out.index,                          
            ch_gtf_for_index.map{it[1]},                     
            GFFREAD.out.gffread_fasta.map{it[1]},            
            false,                                           
            []                                               
        )
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]})
    } else {
        ch_bam_hisat2 = Channel.empty()
        ch_bam_star   = Channel.empty()
        ch_bam_dedup  = Channel.empty()

        if (params.mode == 'hisat2') {
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

            ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.collect{it[1]})

            ch_bam_hisat2 = HISAT2_ALIGN.out.bam
        }
        if (params.mode == 'star') {
            // STAR
            println "Using STAR for alignment"

            if (params.trimmer != "none") {
                ch_samplesheet
                    .map { meta, reads ->
                        def invalid = reads.any { it.toString().endsWith('.gz') }
                        if (invalid) {
                            throw new IllegalArgumentException("❌ gzipped (.gz) fastq files are not yet supported when using a trimmer with star alignment: ${reads}")
                        }
                        tuple(meta, reads)
                }
            }

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
                "n.a."
            )
            ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

            ch_bam_star = STAR_ALIGN.out.bam
        }

        ch_bam = ch_bam_hisat2.mix(ch_bam_star)

        if (params.mark_duplicates) {
            SAMTOOLS_FAIDX (
                ch_fasta_for_index,
                [[],[]],
                false
            )

            SAMTOOLS_SORT(
                ch_bam,
                ch_fasta_for_index,
                "bai"
            )
            ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

            PICARD_MARKDUPLICATES (
                SAMTOOLS_SORT.out.bam,
                ch_fasta_for_index,
                SAMTOOLS_FAIDX.out.fai
            )

            ch_bam_dedup = PICARD_MARKDUPLICATES.out.bam

            ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})
        }

        ch_final_bam = params.mark_duplicates ? ch_bam_dedup : ch_bam

        ch_featurecounts_input = ch_final_bam.map { meta, bam_file ->
            tuple(meta, bam_file, file(ch_gtf))
        }
        SUBREAD_FEATURECOUNTS(ch_featurecounts_input)

        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})
    }

    ch_bam = ch_bam_hisat2.mix(ch_bam_star)

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

    emit:
    multiqc_report      = MULTIQC.out.report.toList()
    versions            = ch_versions                 
    quantification      = params.mode == 'salmon' ? SALMON_QUANT.out.results : SUBREAD_FEATURECOUNTS.out.counts 
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
