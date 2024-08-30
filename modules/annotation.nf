process SNPEFF {

publishDir "${params.output}", mode: 'copy'

input:
path(trimmed_vcf)

output:
path("${params.batch_name}_snpeff.vcf"), emit: annotated_vcf

script:
"""
#snpeff build -genbank -v NC_045512.2
snpeff ann -formatEff -classic -noStats -noLog -quiet -no-upstream -no-downstream NC_045512.2 -v ${trimmed_vcf} > "${params.batch_name}_snpeff.vcf"

"""

}