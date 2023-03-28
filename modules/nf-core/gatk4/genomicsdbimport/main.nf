process GATK4_GENOMICSDBIMPORT {
//    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    path(gvcf_map)
    path fasta
    path fai
    path dict
    path list

    output:
    path("my_database")        , optional:true, emit: genomicsdb
    path("final.vcf.gz") , emit: final_vcfgz
    path("output_filtered.vcf"), emit: filtered_vcfgz
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
//    def args = task.ext.args   ?: ''
//    prefix   = task.ext.prefix ?: "${meta.id}"

    // settings for running default create gendb mode
//    input_command = input_map ? "--sample-name-map ${vcf[0]}" : vcf.collect(){"--variant $it"}.join(' ')

//    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
//    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
//    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
//    if (run_intlist) {
//        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
//        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
//    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
//    if (run_updatewspace) {
//        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
//        interval_command = ''
//        updated_db = "${wspace}"
//    }

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GenomicsDBImport] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" GenomicsDBImport \\
        --sample-name-map $gvcf_map \\
        --genomicsdb-workspace-path my_database \\
        --tmp-dir . \\
        -L $list

    gatk --java-options "-Xmx${avail_mem}g" GenotypeGVCFs \\
        -R $fasta \\
        -V gendb://my_database \\
	-L $list \\
        -O final.vcf.gz

    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
      -V final.vcf.gz \\
      --filter-expression "QUAL < 10.4139" \\
      --filter-name "DRAGENHardQUAL" \\
      -L $list \\
      -O output_filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    genomicsdb_command = "--genomicsdb-workspace-path ${prefix}"
    interval_command = interval_file ? "--intervals ${interval_file}" : "--intervals ${interval_value}"
    updated_db = ""

    // settings changed for running get intervals list mode if run_intlist is true
    if (run_intlist) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = "--output-interval-list-to-file ${prefix}.interval_list"
    }

    // settings changed for running update gendb mode. input_command same as default, update_db forces module to emit the updated gendb
    if (run_updatewspace) {
        genomicsdb_command = "--genomicsdb-update-workspace-path ${wspace}"
        interval_command = ''
        updated_db = "${wspace}"
    }

    def stub_genomicsdb = genomicsdb_command == "--genomicsdb-workspace-path ${prefix}" ? "touch ${prefix}" : ""
    def stub_interval   = interval_command == "--output-interval-list-to-file ${prefix}.interval_list" ? "touch ${prefix}.interval_list" : ""
    def stub_update     = updated_db != "" ? "touch ${wspace}" : ""

    """
    ${stub_genomicsdb}
    ${stub_interval}
    ${stub_update}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
