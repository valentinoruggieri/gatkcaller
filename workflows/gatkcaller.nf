/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGatkcaller.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.list) { ch_list = file(params.list) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
fasta                      = params.fasta ? Channel.fromPath(params.fasta).collect() : Channel.empty()
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome'
include { MAP_GVCF } from '../modules/local/map_gvcf.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { DRAGMAP_ALIGN               } from '../modules/nf-core/dragmap/align/main'
include { GATK4_MARKDUPLICATES        } from '../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_CALIBRATEDRAGSTRMODEL } from '../modules/nf-core/gatk4/calibratedragstrmodel/main'
include { GATK4_HAPLOTYPECALLER       } from '../modules/nf-core/gatk4/haplotypecaller/main_par'
include { GATK4_GENOMICSDBIMPORT      } from '../modules/nf-core/gatk4/genomicsdbimport/main'
// include { GATK4_GATHER_VCF            } from '../modules/local/gathervcf/main'
// include { VCFLIST                     } from '../modules/local/vcflist/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow GATKCALLER {

ch_versions = Channel.empty()


INPUT_CHECK (                       // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    ch_input
)
ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


FASTP (
    INPUT_CHECK.out.reads          // MODULE: Run FastP
)


FASTQC (                           // MODULE: Run FastQC
    FASTP.out.reads
)
ch_versions = ch_versions.mix(FASTQC.out.versions.first())

CUSTOM_DUMPSOFTWAREVERSIONS (       // MODULE: Add description here
    ch_versions.unique().collectFile(name: 'collated_versions.yml')
)


PREPARE_GENOME (                  // SUBWORKFLOW: Prepare genome, samtools index, gatk dictionary and dragmap hash table
    fasta
)


ch_bam_mapped = Channel.empty()

DRAGMAP_ALIGN (                   // MODULE: Run DRAGMAP
    FASTP.out.reads,
    PREPARE_GENOME.out.hashtable,
    'sort'
)


GATK4_MARKDUPLICATES (           // MODULE: Run Markduplicates
    DRAGMAP_ALIGN.out.bam,
    fasta,
    PREPARE_GENOME.out.fasta_fai
)


GATK4_CALIBRATEDRAGSTRMODEL (    // MODULE: Dragon Calibrate STR Model
    GATK4_MARKDUPLICATES.out.bam,
    fasta,
    PREPARE_GENOME.out.fasta_fai,
    PREPARE_GENOME.out.dict,
    PREPARE_GENOME.out.str_file
)

mark = GATK4_MARKDUPLICATES.out.bam.view()

bammodel = Channel.empty().mix(GATK4_MARKDUPLICATES.out.bam, GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model )
        .groupTuple()

bam = Channel.empty().mix(GATK4_MARKDUPLICATES.out.bam, GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model )
        .groupTuple()
        .map { meta, bammodel ->  tuple(meta, bammodel[0])}

model = Channel.empty().mix(GATK4_MARKDUPLICATES.out.bam, GATK4_CALIBRATEDRAGSTRMODEL.out.dragstr_model )
        .groupTuple()
        .map { meta, bammodel -> tuple(meta, bammodel[1])}

//mark = GATK4_MARKDUPLICATES.out.bam.view()
//bam_view = bam.view()
//model_view = model.view()


GATK4_HAPLOTYPECALLER (                // MODULE: HaplotypeCaller
    bam,
    model,
    fasta,
    PREPARE_GENOME.out.fasta_fai,
    PREPARE_GENOME.out.dict,
    [],
    [],
    ch_list
)


// collect map file (sample sample.g.vcf.gz) for GenomicsDBImport
genotype_gvcf_to_call = Channel.empty().mix(GATK4_HAPLOTYPECALLER.out.gvcf_map)
                        .collect()

MAP_GVCF (
        genotype_gvcf_to_call
)

GATK4_GENOMICSDBIMPORT (     // MODULE: GATK4_GENOMICSDBIMPORT.
    MAP_GVCF.out.map_gvcf,
    fasta,
    PREPARE_GENOME.out.fasta_fai,
    PREPARE_GENOME.out.dict,
    ch_list
    )

/*
ch_test = GATK4_HAPLOTYPECALLER.out.vcf.view()

VCFLIST (
    GATK4_HAPLOTYPECALLER.out.vcf
)

ch_test2 = VCFLIST.out.listvcf.view()

GATK4_GATHER_VCF (
    VCFLIST.out.listvcf
)
*/


//
// MODULE: MultiQC
//
    workflow_summary    = WorkflowGatkcaller.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowGatkcaller.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
