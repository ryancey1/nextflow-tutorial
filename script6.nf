// Step 6 - MultiQC report

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"

// exercise 1.1
// Modify the `script1.nf` adding a fourth parameter named `outdir` and set it to a default path
// that will be used as the pipeline output directory.
params.outdir = "$baseDir/out_dir"

// println "reads: $params.reads"

// exercise 1.2
// Modify the script1.nf to print all the pipeline parameters using log.info instead of the println command and a multiline string statement.
// Tip: see an example here: https://github.com/nextflow-io/rnaseq-nf/blob/42974a2/main.nf#L34-L40
log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcript}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """.stripIndent()

/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome from params.transcript
     
    output:
    path 'index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

// Exercise 3.1
// Use the set operator in place of = assignment to define the read_pairs_ch channel.
// read_pairs_ch = Channel .fromFilePairs(params.reads)

// Exercise 3.2
// Use the checkIfExists for the fromFilePairs method to make sure it returns some file pairs

// Exercise 5.1
// Modify the creation of the read_pairs_ch channel by
// using an "into" operator in place of "set"
Channel 
        .fromFilePairs(params.reads, checkIfExists: true)
        // .set { read_pairs_ch }
        .into { read_pairs_ch1; read_pairs_ch2 } // exercise 5.1

/*
 * Run Salmon to perform the quantification of expression using
 * the index and the matched read files
 */
process quantification {

    // Exercise 4.1
    // Add a publishDir directive to store process results into
    // directory of your choice
    publishDir "${params.outdir}/quant"

    input:
    // index_ch output from index process
    path index from index_ch
    // from channel declaration above
    tuple val(pair_id), path(reads) from read_pairs_ch1
 
    output:
    path(pair_id) into quant_ch
 
    script:
    // reads is indexable (reads[0] = first read, reads[1] = second read)
    """
    salmon quant --threads ${task.cpus} \
                 --libType=U \
                 -i ${index} \
                 -1 ${reads[0]} \
                 -2 ${reads[1]} \
                 -o ${pair_id}
    """.stripIndent()
}

/*
 * Run fastQC to check quality of reads files
 */
process fastqc {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads) from read_pairs_ch2

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch
    
    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  


 
/*
 * Create a report using multiQC for the quantification
 * and fastqc processes
 */
process multiqc {
    publishDir params.outdir, mode:'copy'

    // In this script note the use of the mix and collect operators 
    // chained together to get all the outputs of the quantification 
    // and fastqc process as a single input.
    // collect gathers all outputs from fastqc
    // mix combines all collected outputs into one channel       
    input:
    path('*') from quant_ch.mix(fastqc_ch).collect()
    
    output:
    path('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 

// Execute the script with the following command:
// $ nextflow run script6.nf -resume --reads 'data/ggal/*_{1,2}.fq' 

// difference between mix and collect operators
Channel 
    .from( 'a', 'b' )
    .into { c1_1; c1_2 }

Channel 
    .from( 'z' )
    .into { c2_1; c2_2 }

Channel 
    .from( 1, 2, 3 )
    .into { c3_1; c3_2 }

c3_1.mix(c2_1, c1_1).set{ c4 }
println "\nMix only:\n${c4.view()}"

c3_2.mix(c2_2, c1_2).collect().set { c5 }
println "\nMix & collect:\n${c5.view()}"

// mix creates a queue, whereas collect takes a channel and creates a list