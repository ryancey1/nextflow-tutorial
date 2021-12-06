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

println "reads: $params.reads"

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

// Exercise 2.1
// Enable the Docker execution by default adding the below setting in the nextflow.config file.
// docker.enabled = true

// Exercise 2.2
// Print the output of the index_ch channel by using the view operator.
index_ch.view()

// Exercise 2.3 
// Use the command tree work to see how Nextflow organises the process work directory.

// ============= Recap =============
// In this step you have learned:
//     1. How to define a process executing a custom command
//     2. How process inputs are declared
//     3. How process outputs are declared
//     4. How to access the number of available CPUs
//     5. How to print the content of a channel