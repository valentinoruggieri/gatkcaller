

include { SAMPLESHEET_CHECK } from './modules/local/samplesheet_check'

ch_samplesheet = Channel.fromPath(params.input)

workflow INPUT_CHECK {

    SAMPLESHEET_CHECK ( ch_samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }
	.view()
}

workflow {
	INPUT_CHECK()
}
