// Step 1 - define the pipeline parameters

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

// ============= Recap =============
// In this step you have learned:
//     1. How to define parameters in your pipeline script
//     2. How to pass parameters by using the command line
//     3. The use of $var and ${var} variable placeholders
//     4. How to use multiline strings
//     5. How to use log.info to report values