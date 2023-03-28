// params.list = '' // path to list file
process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.3.0.0 conda-forge::parallel bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(model)
    // path(input), path(input_index), path(intervals), 7
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  dbsnp_tbi
    path list

    output:
    tuple val(meta), path("*.g.vcf.gz")       , emit: gvcf
    tuple val(meta), path("*.tbi")          , optional:true, emit: tbi
    tuple val(meta), path("*.realigned.bam"), optional:true, emit: bam
    tuple val(meta), path("*.list")         , emit: vcflist
    tuple val(meta), path("*.gathered.g.vcf.gz")  , emit: gathered_gvcf
    tuple val(meta), path("*.gathered.g.vcf.gz.tbi"), emit: gathered_tbi
    path("*.gvcf.list.map"), emit: gvcf_map
//    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
//    def markedbam_command = "${bammodel[1].map { it ->  (it.first())}}"
//    def model_command = "${bammodel[1].map { it ->  (it.last())}}"
//    def model_command = "${bammodel[2]}"
//    def interval_command = intervals ? "--intervals $intervals" : ""
//    def dragstr_command = dragstr_model ? "--dragstr-params-path !{dragstr_model}" : ""
//    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output !{prefix.replaceAll('.g\\s*!', '')}.realigned.bam" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    parallel -a $list 'myvar={}; gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
        --input ${bam} \\
        --output ${prefix}.\$myvar.g.vcf.gz \\
        --reference ${fasta} \\
        --dragstr-params-path ${model}  \\
        -L \$myvar \\
        -ERC GVCF \\
        --native-pair-hmm-threads 1'

        awk -v pre="${prefix}" '{print pre"."\$1".g.vcf.gz"}' $list > ${prefix}.gathered.list
        
        gatk --java-options "-Xmx${avail_mem}g" GatherVcfs -I ${prefix}.gathered.list -O ${prefix}.gathered.g.vcf.gz
        tabix -p vcf ${prefix}.gathered.g.vcf.gz
        realpath ${prefix}.gathered.g.vcf.gz | awk -v pre="${prefix}" '{print pre"\t"\$1}' > ${prefix}.gvcf.list.map

    """
}

/*
    awk -v pre="${prefix}" '{print pre"."\$1".vcf.gz"}' /mnt/c/Users/Public/Dropbox/Pipelines/Workflow/nextflow/nf-core-gatkcaller/list > ${prefix}.gathered.list
//        realpath ${prefix}.*.vcf.gz > ${prefix}.vcf.list
    gatk --java-options "-Xmx${avail_mem}g" GatherVcfs
        -I ${prefix}.gathered.list
        -O ${prefix}.gathered.vcf.gz
        --CREATE_INDEX
        /mnt/c/Users/Public/Dropbox/Pipelines/Workflow/nextflow/nf-core-gatkcaller/list
*/