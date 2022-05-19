#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*  Take the first N reads of a pair of FASTQ files and output truncated FASTQ files
 *  This provides a simple way to get a small testset of reads
 *  based on:
 *  https://stackoverflow.com/questions/69010242/nextflow-how-should-i-truncate-a-fastq-file-at-line-x-process-fails-with-error
 *  usgage example:
 *    nextflow run readpickn/main.nf -with-docker [nextflow-image] --datadir [path or s3 folder uri] --outdir [path or s3 folder uri]
 */ 

params.max_records = 10000
params.max_lines = 4*params.max_records
params.datadir = '.'
params.outdir = 'results'

process TRUNCATE_FASTQ {

    container 'quay.io/biocontainers/pysam:0.16.0.1--py37hc334e0b_1'
    publishDir "$params.outdir"    

    input:
    tuple val(sample), path(reads)

    output:
    path "*.fastq.gz"

    script:
    def (fq1, fq2) = reads
    println fq1
    println fq2      

    """
    #!/usr/bin/env python
    import gzip
    import pysam

    def truncate(fn, num):
        with pysam.FastxFile(fn) as ifh:
            with gzip.open(f'truncated_{fn}', 'wt') as ofh:
                for idx, entry in enumerate(ifh):
                    if idx >= num:
                        break
                    print(str(entry), file=ofh)

    truncate('${fq1}', ${params.max_records})
    truncate('${fq2}', ${params.max_records})
    truncate('${fq2}', ${params.max_records})
    """
}      

process TRUNCATE_FASTQ_SIMPLE {
    publishDir "$params.outdir"

    input:
    tuple val(sample), path(reads)

    output:
    path "*.fastq.gz"

    script:
    def (fq1, fq2) = reads

    """
    head -n ${params.max_lines} < <(gzcat "${fq1}") | gzip > "truncated_${fq1}"
    head -n ${params.max_lines} < <(gzcat "${fq2}") | gzip > "truncated_${fq2}"
    """
}

process printInputs {
    input:
    tuple val(sample), path(reads)

    script:
    def (fq1, fq2) = reads      

    exec:
    println "first input file: $fq1"
    println "second input file: $fq2"
}

workflow {

    globstr = params.datadir + '/*{1,2}.fastq.gz'
    println "glob string: " + globstr
    inputs = Channel.fromFilePairs(globstr)
    printInputs( inputs )

    TRUNCATE_FASTQ( inputs )
//    TRUNCATE_FASTQ_SIMPLE( inputs )
}
