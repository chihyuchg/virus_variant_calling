process TRIMMOMATIC {
    
publishDir "${params.output}", mode: 'copy'

input:
tuple val(sample_id), path(fq_r1), path(fq_r2)

output:
tuple val(sample_id), path("${sample_id}_trimmed_r1.fq"), path("${sample_id}_trimmed_r2.fq"), emit: trimmed_fqs
tuple val(sample_id), path("trim_log.txt"), emit: trim_log

script:
"""
trimmomatic PE \
-phred33 \
${fq_r1} ${fq_r2} \
${sample_id}_trimmed_r1.fq ${sample_id}_trimmed_unpaired_r1.fq \
${sample_id}_trimmed_r2.fq ${sample_id}_trimmed_unpaired_r2.fq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
2> trim_log.txt

"""
}