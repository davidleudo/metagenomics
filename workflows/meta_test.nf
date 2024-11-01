/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW                                  } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                               } from '../modules/nf-core/multiqc/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_reads = Channel.of([
                        [ id:'meta_test', single_end:false ], // meta map
                        [ 
                            file( 'https://drive.google.com/file/d/1xic_nvxwWRrRkxk-1YCpImW_Z49i8CuI/view?usp=drive_link', checkIfExists: true ), 
                            file( 'https://drive.google.com/file/d/17B5DACcQj-UGsAa2_jPYDg3YLACtJrE_/view?usp=drive_link', checkIfExists: true )
                        ]
                    ])

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
}

workflow {
    META_TEST()
}