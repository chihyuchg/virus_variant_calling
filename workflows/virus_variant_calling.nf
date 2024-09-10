#!/usr/bin/env nextflow

/*
=========================================================================
    IMPORT MODULES/SUBWORKFLOWS
=========================================================================
*/

include {TRIMMOMATIC        } from '../modules/trim_fq'
include {HISAT2             } from '../modules/alignment'
include {HAPLOTYPECALLER    } from '../modules/gatk'
include {GENOTYPEGVCFS        } from '../modules/gatk'
include {VARIANTFILTER      } from '../modules/gatk'  
include {LEFTALIGNANDTRIMVARIANTS} from '../modules/gatk'
include {GATKREF           } from '../modules/gatk' 
include {SNPEFF             } from '../modules/annotation'

/*
=========================================================================
    WORKFLOW DEFINITIONS
=========================================================================
*/

workflow VIRUS_VARIANT_CALLING {
    take:
        sample_ch
        ref_fa


    main:
        TRIMMOMATIC (sample_ch)
        HISAT2 (ref_fa, TRIMMOMATIC.out.trimmed_fqs)
        HAPLOTYPECALLER(ref_fa, HISAT2.out.bam)

    emit:
        gvcfs = HAPLOTYPECALLER.out.gvcf
        gvcfs_idx = HAPLOTYPECALLER.out.gvcf_idx

}

workflow COMBINE_VCFS {
    take:
        gvcf_ch
        gvcf_idx_ch
        ref_fa

    main:
        GATKREF (ref_fa, gvcf_ch)
        GENOTYPEGVCFS (GATKREF.out.refs, gvcf_ch, gvcf_idx_ch)
        VARIANTFILTER (GATKREF.out.refs, GENOTYPEGVCFS.out.combined_vcf)
        LEFTALIGNANDTRIMVARIANTS (GATKREF.out.refs, VARIANTFILTER.out.filtered_vcf)
        SNPEFF (LEFTALIGNANDTRIMVARIANTS.out.trimmed_vcf)


    emit:
        vcf = SNPEFF.out.annotated_vcf

}


workflow.onError {
    println "Opps... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}