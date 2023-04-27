params.read_triplets = "$projectDir/test_data/reads/reads_{1,2,3}*.gz"
params.outdir = "results"

process UMI_TRANSFER {
    input:
    tuple val(sample_id), path(read_triplets)

    output:
    tuple val(sample_id), path('*_with_umi.fastq.gz')

    publishDir params.outdir
    script:
    /*
    First decompress and put 12 fastq rows (1 entry from each file) in one row as a tsv
    awk handles the UMI insertion into the header row of the read 1 and 3.
    Details to read the awk command:
    UMI should be placed right before the first space of column \$1 (and \$3 respectively).
    Column \$5 will be the UMI sequence (second row of file 2)
    \${sample_id}_1_with_umi.fastq.gz will be built up by column \$1, \$4, a '+' and \$10.
    \${sample_id}_2_with_umi.fastq.gz will be built up by column \$3, \$6, a '+' and \$12.
    */
    """
        paste <(gzip -dc ${read_triplets[0]}) <(gzip -dc ${read_triplets[1]}) <(gzip -dc ${read_triplets[2]}) | paste - - - - | \
        awk -F "\\t" 'BEGIN{OFS="\\n"}; \
            match(\$1,/@[A-Za-z0-9:_]+/) {print substr( \$1, RSTART, RLENGTH )":"\$5substr( \$1, RSTART+RLENGTH ),\$4,"+",\$10 | "gzip > ${sample_id}_1_with_umi.fastq.gz"}; \
            match(\$3,/@[A-Za-z0-9:_]+/) {print substr( \$3, RSTART, RLENGTH )":"\$5substr( \$3, RSTART+RLENGTH ),\$6,"+",\$12 | "gzip > ${sample_id}_2_with_umi.fastq.gz"}'
    """
}


pattern = ~/(\w+)_R*[1-3]\.(fastq|fq).gz$/

workflow {
    Channel.fromPath(params.read_triplets).view()
    Channel
        .fromFilePairs(params.read_triplets, size: 3) { file -> (file.name =~ pattern)[0][1] }
        .set{ read_triplets_ch }
    read_triplets_ch.view()
    UMI_TRANSFER(read_triplets_ch)
}