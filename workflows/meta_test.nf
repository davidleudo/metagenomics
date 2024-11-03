/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW                                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM                                 } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                                           } from '../modules/nf-core/trimmomatic/main'
include { MULTIQC                                               } from '../modules/nf-core/multiqc/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
ch_reads = Channel.of([
                        [ id:'meta_test', single_end:false ], // meta map
                        [ 
                            file( '/workspace/metagenomics/assets/data/muestra09.R1-003.fastq.gz', checkIfExists: true ), 
                            file( '/workspace/metagenomics/assets/data/muestra09.R2-004.fastq.gz', checkIfExists: true )
                        ]
                    ])
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow META_TEST {

    take:
    ch_reads

    main:
    ch_versions = Channel.empty()

    FASTQC_RAW(
        ch_reads
    )
    
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    TRIMMOMATIC(
        ch_reads
    )
    ch_versions= ch_versions.mix(TRIMMOMATIC.out.versions.first())

    
}

/*
workflow {
    META_TEST(ch_reads)
}
*/