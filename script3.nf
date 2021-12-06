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

// Exercise 3.1
// Use the set operator in place of = assignment to define the read_pairs_ch channel.
// read_pairs_ch = Channel .fromFilePairs(params.reads)

// Exercise 3.2
// Use the checkIfExists for the fromFilePairs method to make sure it returns some file pairs
Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

// Edit the script script3.nf and add the following statement as the last line:
// read_pairs_ch.view()
read_pairs_ch.view()

// specify different read files by using a glob pattern
// nextflow run script3.nf --reads 'data/ggal/*_{1,2}.fq'
// output:  [liver, [/home/ryancey/git/nextflow-tutorial/data/ggal/liver_1.fq, /home/ryancey/git/nextflow-tutorial/data/ggal/liver_2.fq]]
//          [lung, [/home/ryancey/git/nextflow-tutorial/data/ggal/lung_1.fq, /home/ryancey/git/nextflow-tutorial/data/ggal/lung_2.fq]]
//          [gut, [/home/ryancey/git/nextflow-tutorial/data/ggal/gut_1.fq, /home/ryancey/git/nextflow-tutorial/data/ggal/gut_2.fq]]

// ============= Recap =============
// In this step you have learned:
//     1. How to use fromFilePairs to handle read pair files.
//     2. How to use the set operator to define a new channel variable.
//     3. How to use the checkIfExists option.