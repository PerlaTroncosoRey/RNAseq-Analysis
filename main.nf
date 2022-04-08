#!/usr/bin/env nextflow

/*
 *	RNAseq pipeline implemented with Nextflow
 *	Author:
 *	Perla Rey <Perla.Rey at quadram.ac.uk>
*/


timestamp='20220330'

/*
 * Defines some parameters in order to specify the reference genome
 * and read pairs by using the command line options
 */
params.version='0.0.1'
params.reads = "$projectDir/data/nat-prot-example/chrX_data/samples/*_{1,2}.fastq.gz"
params.hisat2indexes = "$projectDir/data/nat-prot-example/chrX_data/indexes/chrX_tran"
params.genes = "$projectDir/data/nat-prot-example/chrX_data/genes/chrX.gtf"
// params.annot = "$projectDir/data/nat-prot-example/chrX_data/ensembl/.bed.gff"
// params.genome = "$projectDir/data/nat-prot-example/chrX_data/.fa"
params.outdir = 'results'


/*
 * Prints version when asked for
 */

if (params.version) {
	System.out.println("")
	System.out.println("RNA-seq Pipeline with Nextflow - Version $params.version ($timestamp)")
	//exit 1
}


/*
 * Create the reads_pair_ch channel containing three elements:
 * SampleID, forward read (R1), reverse read (R2)
 */

Channel
     .fromFilePairs( params.reads )
     .ifEmpty { error "Cannot find any reads matching: ${params.reads}"}
     .set { read_pairs_ch }



process qcFastp{
  input:
  tuple val(pair_id), path(reads) from read_pairs_ch

  output:
  set pair_id, "fq.gz" into fqFastp_ch
  set pair_id, "html" into htmlFastp_ch
  set pair_id, "json" into jsonFastp_ch


  """
  fastp --version
  fastp -i ${pair_id}_1.fastq.gz -o ${pair_id}_1.fq.gz \
    -I ${pair_id}_2.fastq.gz -O ${pair_id}_2.fq.gz \
    --detect_adapter_for_pe -q 30 -l 60 -w ${task.cpus} \
		-h ${pair_id}.fastp.html -j ${pair_id}.fastp.json
  """
}
