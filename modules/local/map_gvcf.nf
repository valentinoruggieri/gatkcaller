process MAP_GVCF {
//    tag "$meta.id"
    label 'process_medium'

//    conda "bioconda::gatk4=4.3.0.0 conda-forge::parallel bioconda::tabix=1.11"
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
//        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    path(gvcf)

    output:
    path("map_gvcf.map"),      emit: map_gvcf

    script:
"""
cat $gvcf > map_gvcf.map

"""
}

// | awk '{print \$1"\t"\$1}'
