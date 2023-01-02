#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """
    if [[ "${params.run_type}" == "r2d2" ]] || [[ "${params.run_type}" == "raven" ]] ; 
      then
        cd ${params.image_folder}
        if [[ ! -f deepvariant-1.4.0.sif ]] ;
          then
            singularity pull deepvariant-1.4.0.sif docker://google/deepvariant:1.4.0
        fi
    fi
    if [[ "${params.run_type}" == "local" ]] ; 
      then
        docker pull google/deepvariant:1.4.0
    fi
    """

}

process deepvariant {
  tag "${f}"
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    path f
  
  script:
    """
    mkdir -p /workdir/deepvariant_output
    
    if [ "${params.exomebed}" != ""  ] ; then
    singularity run --cleanenv \
                    --no-home \
                    -B ${params.tmp}:/workdir/tmp \
                    ${params.image_folder}/deepvariant-1.4.0.sif \
                    /opt/deepvariant/bin/run_deepvariant \
                    --model_type=${params.model} \
                    --ref=${params.genome}/${params.organism}.${params.release}.fa \
                    --reads=${params.deepvariant_raw_data}/${f} \
                    --regions=${params.exomebed} \
                    --output_vcf=/workdir/deepvariant_output/${f%${.sorted.bam}}.vcf.gz \
                    --output_gvcf=/workdir/deepvariant_output/${f%${.sorted.bam}}.g.vcf.gz \
                    --sample_name ${f%${.sorted.bam}} \
                    --num_shards=18
    else
    singularity run --cleanenv \
                    --no-home \
                    -B ${params.tmp}:/workdir/tmp \
                    ${params.image_folder}/deepvariant-1.4.0.sif \
                    /opt/deepvariant/bin/run_deepvariant \
                    --model_type=${params.model} \
                    --ref=${params.genome}/${params.organism}.${params.release}.fa \
                    --reads=${params.deepvariant_raw_data}/${f} \
                    --output_vcf=/workdir/deepvariant_output/${f%${.sorted.bam}}.vcf.gz \
                    --output_gvcf=/workdir/deepvariant_output/${f%${.sorted.bam}}.g.vcf.gz \
                    --sample_name ${f%${.sorted.bam}} \
                    --num_shards=18
    fi
    """
}

workflow images {
  main:
    get_images()
}


workflow {
    data = channel.fromPath( "${params.deepvariant_raw_data}/*.sorted.bam" )
    // the following command makes sure it is only run when an output file is missing
    // data = data.filter{ ! file("$it".replaceAll(/.fastq.gz/, "_fastqc.html").replace("${params.fastqc_raw_data}", "${params.project_folder}/fastqc_output/") ).exists() }
    deepvariant( data )
}