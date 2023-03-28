//
// Mark duplicate and calibrate dragmap model
//


include { GATK4_MARKDUPLICATES        } from '../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_CALIBRATEDRAGSTRMODEL } from '../modules/nf-core/gatk4/calibratedragstrmodel/main'


workflow MARK_CALIBRATE {
    take:
        bam
        fasta
        fai
        dict
        str_file                  // channel: [mandatory] fasta

    main:

    ch_versions = Channel.empty()


    GATK4_MARKDUPLICATES (           // MODULE: Run Markduplicates
    bam,
    fasta,
    fai
)


GATK4_CALIBRATEDRAGSTRMODEL (    // MODULE: Dragon Calibrate STR Model
    GATK4_MARKDUPLICATES.out.bam,
    fasta,
    fai,
    dict,
    str_file
)

    emit:
