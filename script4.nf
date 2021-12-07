// Step 4 - Perform expression quantification

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"
params.transcript = "$baseDir/data/ggal/transcriptome.fa"
params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcript}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

 
/* 
 * define the `index` process that create a binary index 
 * given the transcriptome file
 */
process index {
    
    input:
    path transcriptome from params.transcript
     
    // note the index_ch output 
    output:
    path 'index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


Channel 
    .fromFilePairs( params.reads, checkIfExists:true )
    .set { read_pairs_ch } 

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
    tuple val(pair_id), path(reads) from read_pairs_ch
 
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

// try to execute it with more read files as shown below:
// 
// $ nextflow run script4.nf -resume --reads 'data/ggal/*_{1,2}.fq'
// 
// the quantification process is executed more than one time.

// ============= Recap =============
// In this step you have learned:
//    1. How to connect two processes by using the channel declarations
//    2. How to resume the script execution skipping already already computed steps
//    3. How to use the publishDir to store a process results in a path of your choice