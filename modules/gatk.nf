process GATKREF {

publishDir "${params.output}", mode: 'copy'

input:
path ref_fa
path gvcfs

output:
tuple env(ref_name), path(ref_fa), path("*.fai"), path("*.dict"), emit: refs


script:
"""
ref=${ref_fa}
ref_name=\${ref%.*}

samtools faidx ${ref_fa}
gatk CreateSequenceDictionary -R ${ref_fa}

"""


}


process HAPLOTYPECALLER {

publishDir "${params.output}", mode: 'copy'

input:
path ref_fa
tuple val(sample_id), path(input_bam), path(input_bai)

output:
path("${sample_id}.g.vcf"), emit: gvcf
path("${sample_id}.g.vcf.idx"), emit: gvcf_idx



script:
"""
samtools faidx ${ref_fa}
gatk CreateSequenceDictionary -R ${ref_fa}

gatk HaplotypeCaller \
-R ${ref_fa} \
-I ${input_bam} \
-O "${sample_id}.g.vcf" \
--minimum-mapping-quality 10 \
--ploidy 2 -ERC BP_RESOLUTION

"""

}


process COMBINEVCFS {

publishDir "${params.output}", mode: 'copy'

input:
tuple val(ref_name), path(ref_fa), path("${ref_name}.fna.fai"), path("${ref_name}.dict")
path(gvcfs)
path(gvcfs_idx)

output:
path("${params.batch_name}_combined.vcf"), emit: combined_vcf

script:
def input_list = gvcfs.collect{"--variant $it"}.join(' ')

"""

gatk CombineGVCFs \
-R ${ref_fa} \
${input_list} \
-O "${params.batch_name}_combined.vcf"

"""

}


process VARIANTFILTER {

publishDir "${params.output}", mode: 'copy'

input:
tuple val(ref_name), path(ref_fa), path("${ref_name}.fna.fai"), path("${ref_name}.dict")
path(combined_vcf)

output:
path("${params.batch_name}_filtered.vcf"), emit: filtered_vcf


script:
"""

gatk VariantFiltration \
-R ${ref_fa} \
-V ${combined_vcf} \
-O "${params.batch_name}_filtered.vcf" \
--filter-name "lowQUAL100" \
--filter-expression "QUAL < 100" \
--filter-name "lowQD2.0" \
--filter-expression "QD < 2.0" \
--filter-name "lowReadPosRankSum4.0" \
--filter-expression "ReadPosRankSum < -4.0"

"""

}


process LEFTALIGNANDTRIMVARIANTS {

input:
tuple val(ref_name), path(ref_fa), path("${ref_name}.fna.fai"), path("${ref_name}.dict")
path(filtered_vcf)

output:
path("${params.batch_name}_trimmed_variants.vcf"), emit: trimmed_vcf

script:
"""

gatk LeftAlignAndTrimVariants \
--split-multi-allelics \
-R ${ref_fa} \
-V ${filtered_vcf} \
-O "${params.batch_name}_trimmed_variants.vcf"

"""

}
