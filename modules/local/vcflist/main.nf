process VCFLIST {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.3.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(gzvcf)

    output:
    path ("*.gzvcf.list")                                        , emit: listvcf


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
//    def avail_mem = 3
//    if (!task.memory) {
//        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
//    } else {
//        avail_mem = task.memory.giga
//    }

"""
    echo $gzvcf | tr ' ' '\\n' > ${prefix}.gzvcf.list
    cat ${prefix}.gzvcf.list

"""
}

