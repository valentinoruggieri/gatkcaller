//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments
// Samtools index fasta | gatk dict | dragmap hash table

include { DRAGMAP_HASHTABLE                      } from '../../modules/nf-core/dragmap/hashtable/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                         } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_COMPOSESTRTABLEFILE              } from '../../modules/nf-core/gatk4/composestrtablefile/main'

workflow PREPARE_GENOME {
    take:
        fasta                   // channel: [mandatory] fasta

    main:

    ch_versions = Channel.empty()

    DRAGMAP_HASHTABLE(fasta.map{ it -> [[id:it[0].baseName], it] }) // If aligner is dragmap
    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    SAMTOOLS_FAIDX(fasta.map{ it -> [[id:it[0].baseName], it] })
    GATK4_COMPOSESTRTABLEFILE(fasta, SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }, GATK4_CREATESEQUENCEDICTIONARY.out.dict)

    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [file1, file2] becomes [[meta1, file1],[meta2, file2]]
    // outputs are collected to maintain a single channel for relevant TBI files
//    TABIX_DBSNP(dbsnp.flatten().map{ it -> [[id:it.baseName], it] })
//    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [[id:it.baseName], it] })
//    TABIX_KNOWN_SNPS( known_snps.flatten().map{ it -> [[id:it.baseName], it] } )
//    TABIX_KNOWN_INDELS( known_indels.flatten().map{ it -> [[id:it.baseName], it] } )
//    TABIX_PON(pon.flatten().map{ it -> [[id:it.baseName], it] })


//    chr_files = chr_dir
//    if (params.chr_dir && params.chr_dir.endsWith('tar.gz')) {
//        UNTAR_CHR_DIR(chr_dir.map{ it -> [[id:it[0].baseName], it] })
//        chr_files = UNTAR_CHR_DIR.out.untar.map{ it[1] }
//        ch_versions = ch_versions.mix(UNTAR_CHR_DIR.out.versions)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_versions = ch_versions.mix(DRAGMAP_HASHTABLE.out.versions)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    emit:
        hashtable                        = DRAGMAP_HASHTABLE.out.hashmap.map{ meta, index -> [index] }.collect() // path: dragmap/*
        dict                             = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta_fai                        = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        str_file                         = GATK4_COMPOSESTRTABLEFILE.out.str_table                               // path: str_table.zip file
//        chr_files                        = chr_files

        versions                         = ch_versions                                                         // channel: [ versions.yml ]
}
