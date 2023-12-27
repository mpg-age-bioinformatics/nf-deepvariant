#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ;

      then

        cd ${params.image_folder}

        if [[ ! -f deepvariant-1.4.0.sif ]] ;
          then
            singularity pull deepvariant-1.4.0.sif docker://google/deepvariant:1.4.0
        fi

        if [[ ! -f rnaseq.python-3.8-2.sif ]] ;
          then
            singularity pull rnaseq.python-3.8-2.sif docker://index.docker.io/mpgagebioinformatics/rnaseq.python:3.8-2
        fi

    fi

    if [[ "${params.containers}" == "docker" ]] ;

      then

        docker pull google/deepvariant:1.4.0
        ocker pull mpgagebioinformatics/rnaseq.python:3.8-2

    fi
    """

}

process ucsc_to_ensembl {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val bed

  output:
    val bed

  when:
    ( "${bed}" != "none")

  script:
    """
    sed -i 's/chr//g' ${bed}
    """
}

process deepvariant {
  stageInMode 'symlink'
  stageOutMode 'move'
  
  input:
    tuple val(pair_id), path(bwa)
    val exomebed

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/deepvariant_output/${pair_id}.vcf.gz").exists() ) 
  
  script:
    """
    mkdir -p /workdir/deepvariant_output
    
    if [[ "${exomebed}" != "" ]] ; then

      /opt/deepvariant/bin/run_deepvariant --model_type=${params.model} \
      --ref=${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa \
      --reads=/raw_data/${pair_id}.sorted.bam \
      --regions=${exomebed} \
      --output_vcf=/workdir/deepvariant_output/${pair_id}.vcf.gz \
      --output_gvcf=/workdir/deepvariant_output/${pair_id}.g.vcf.gz \
      --sample_name ${pair_id} \
      --num_shards=${task.cpus}

    elif [[ "${exomebed}" == "" ]] ; then

      /opt/deepvariant/bin/run_deepvariant --model_type=${params.model} \
      --ref=${params.genomes}${params.organism}/${params.release}/${params.organism}.${params.release}.fa \
      --reads=/raw_data/${pair_id}.sorted.bam \
      --output_vcf=/workdir/deepvariant_output/${pair_id}.vcf.gz \
      --output_gvcf=/workdir/deepvariant_output/${pair_id}.g.vcf.gz \
      --sample_name ${pair_id} \
      --num_shards=${task.cpus}
     
    fi
    """
}

process filtering {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id
    tuple val(pair_id), path(bwa)
    path(bed)

  output:
    val pair_id

  when:
    ( ! file("${params.project_folder}/filter/${pair_id}.SNPs.vcf").exists() ) 

  script:
    """
    mkdir -p /workdir/filter
    cd /workdir/deepvariant_output/
    zcat ${pair_id}.vcf.gz | grep PASS > /workdir/filter/${pair_id}.SNPs.vcf
    """
}

process subtractWT {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val pair_id

  script:
    """
    #!/usr/local/bin/python

    import os
    import pandas as pd
    import numpy as np

    filterfolder="/workdir/filter/"
    print(filterfolder)

    sample_sheet="/workdir/sample_sheet.xlsx"
    print(sample_sheet)

    samples = pd.read_excel(sample_sheet, sheet_name='samples')
    print(samples)

    samples['SampleID'] = samples['Group'] + '.Rep_' + samples['Replicate'].astype(str)
    samples

    # dictionary for mapping sample names to sample files
    sample_dict = samples.set_index('Sample').to_dict()['SampleID']
    for index, sample in samples.iterrows():
        print(sample['Sample'])
        # filter snps from mutated samples if background is present
        if not pd.isnull(sample['Background Sample']) and sample_dict[sample['Sample']] != sample_dict[sample['Background Sample']]:
            # read in sample vcf
            tmp= pd.read_csv(filterfolder + sample_dict[sample['Sample']] + '.SNPs.vcf', header=None, sep="\t", comment="#")
            tmp["ref"]=tmp[0].astype(str)+"_"+tmp[1].astype(str)
            # read in corresponding WT vcf,  
            # here we assume that all variants that passed through google's deepVariant filter are not simply sequencing errors, 
            # so no additional filtering is necessary
            wt=pd.read_csv(filterfolder + sample_dict[sample['Background Sample']] + '.SNPs.vcf', header=None, usecols=[0,1], sep="\t", comment="#")
            wt=list(wt[0].astype(str)+"_"+wt[1].astype(str))
            tmp=tmp[~tmp["ref"].isin(wt)]
            tmp=tmp.drop(['ref'], axis=1)
            tmp.to_csv(filterfolder + sample_dict[sample['Sample']] + '.SNPs.nowt.vcf', sep='\t', index=None, header=None)
        else:
            # only write background to file?
            tmp= pd.read_csv(filterfolder + sample_dict[sample['Sample']] + '.SNPs.vcf', header=None, sep="\t", comment="#")
            tmp.to_csv(filterfolder + sample_dict[sample['Sample']] + '.SNPs.nowt.vcf', sep='\t', index=None, header=None)

    """

}

workflow images {
  main:
    get_images()
}

workflow {
  main:
    if ( "${params.exomebed}" != "none" ) {
    exomebed=ucsc_to_ensembl( "${params.exomebed}" )  
    } else {
    exomebed=""
    }
    println "${exomebed}"

    data = channel.fromFilePairs( "${params.deepvariant_raw_data}/*.sorted.bam", size: -1 )
    deepvariant( data,  "${exomebed}")
    filtering( deepvariant.out.collect(), data )
    subtractWT( filtering.out.collect() )
}