nextflow.enable.dsl=2

// -------------------
// Parameters
// -------------------
params.reads     = "data/sample.fastq.gz"
params.reference = "data/reference.fa"
params.outdir    = "results"

// -------------------
// Processes
// -------------------
process QC {
    tag "$reads"

    input:
    path reads

    output:
    path "qc_report.html"

    script:
    """
    fastqc -o ./ $reads
    mv *.html qc_report.html
    """
}

process ALIGN {
    tag "$reads"

    input:
    path reads
    path ref

    output:
    path "aligned.bam"

    script:
    """
    bwa index $ref
    bwa mem $ref $reads | samtools view -Sb - > aligned.bam
    """
}

process VARIANT_CALL {
    tag "$bam"

    input:
    path bam
    path ref

    output:
    path "variants.vcf"

    script:
    """
    samtools sort -o sorted.bam $bam
    samtools index sorted.bam
    bcftools mpileup -f $ref sorted.bam | bcftools call -mv -Ov -o variants.vcf
    """
}

// -------------------
// Workflow
// -------------------
workflow {
    reads     = file(params.reads)
    reference = file(params.reference)

    qc_report = QC(reads)
    bam_out   = ALIGN(reads, reference)
    variants  = VARIANT_CALL(bam_out, reference)

    // Publish results
    qc_report.view { "QC Report → ${it}" }
    bam_out.view   { "Aligned BAM → ${it}" }
    variants.view  { "Variants VCF → ${it}" }
}
