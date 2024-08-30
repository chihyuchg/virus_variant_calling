process HISAT2 {

publishDir "${params.output}", mode: 'copy'

executor 'local'

input:
path ref_fa
tuple val(sample_id), path(trimmed_fq_r1), path(trimmed_fq_r2)

output:
tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), emit: bam

script:
"""

hisat2-build ${ref_fa} reference
hisat2 -p 2 --no-spliced-alignment --no-unal -x reference -q -1 ${trimmed_fq_r1} -2 ${trimmed_fq_r2} --rg "ID:${sample_id}" --rg "LB:${sample_id}" --rg "SM:${sample_id}" | \
samtools view -Sb -F256 - | \
samtools sort > ${sample_id}_sorted.bam  
samtools index ${sample_id}_sorted.bam

"""

}