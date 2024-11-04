/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_RAW                                  } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM                                 } from '../modules/nf-core/fastqc/main'
include { TRIMMOMATIC                                           } from '../modules/nf-core/trimmomatic/main'
include { MULTIQC                                               } from '../modules/nf-core/multiqc/main'
include { MEGAHIT                                               } from '../modules/nf-core/megahit/main'

//to do Bowtie2 subworkflow
include { BOWTIE2_ASSEMBLY_BUILD                                } from '../modules/local/bowtie2_assembly_build'
include { BOWTIE2_ASSEMBLY_ALIGN                                } from '../modules/local/bowtie2_assembly_align'

//to do Metabat2 subworkflow
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS                  } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { METABAT2_METABAT2                                     } from '../modules/nf-core/metabat2/metabat2/main'

//to do Checkm subworkflow
include { CHECKM_LINEAGEWF                                      } from '../modules/nf-core/checkm/lineagewf/main.nf'

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
    //
    FASTQC_RAW(
        ch_reads
    )
    
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    //
    TRIMMOMATIC(
        ch_reads
    )
    ch_versions= ch_versions.mix(TRIMMOMATIC.out.versions.first())

    //
    FASTQC_TRIM (
        TRIMMOMATIC.out.trimmed_reads
    )
    // incluir canal que separe la salida
    trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    // Canal intermedio para descomponer las rutas en archivos individuales
    separated_reads = trimmed_reads.map { meta, paths ->
        def (reads1, reads2) = paths
        return tuple(meta, reads1, reads2)
    }

    //
    MEGAHIT (separated_reads)

    //
    BOWTIE2_ASSEMBLY_BUILD (MEGAHIT.out.contigs)


    // Combina las salidas para que el siguiente proceso reciba todos los inputs necesarios
    pre_bowtie2align_channel = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index.combine(TRIMMOMATIC.out.trimmed_reads)
    BOWTIE2_ASSEMBLY_ALIGN (pre_bowtie2align_channel)

    //
    pre_metabat2jgi_channel = BOWTIE2_ASSEMBLY_ALIGN.out.mappings.map{ tuple(it[0], it[2], it[3]) }
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS (pre_metabat2jgi_channel)

    // val megahitcontigs depth
    pre_metabat2_channel =  MEGAHIT.out.contigs.combine(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth.map{it[1]})
    //pre_metabat2_channel.view()
    METABAT2_METABAT2 (pre_metabat2_channel)
/*
    CHECKM_LINEAGEWF (
        METABAT2_METABAT2.out.fasta,
        "fa",
        ".")
*/
}
