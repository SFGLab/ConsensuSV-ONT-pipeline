#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and fastq files and can be provided as command line options
 */

params.ref = "/tools/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.outdir = "results"
params.design = "/tools/files.csv"	
params.threads = 4
params.mem = 4

workflow {
    files = Channel.fromPath(params.design).splitCsv()
	ALIGN_PB(files)
	ALIGN_ONT(files)
	Sniffles(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	CuteSV(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Svim(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Dysgu(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Nanovar(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	//NanoSV(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index) execution time too long
	PBSV(ALIGN_PB.out.bam, ALIGN_PB.out.sample, ALIGN_PB.out.index)
}

process ALIGN_PB {
	tag "Create align for Pac Bio Callers"

	publishDir "${params.outdir}/aligment/${fastq.simpleName}", mode: "copy"
	
	input:
	path fastq
	
	output:
	path 'output_pb.bam', emit: bam
	path 'output_pb.bam.bai', emit: index
	val fastq.simpleName, emit: sample
	
	script:
	"""
	pbmm2 align ${params.ref} $fastq output_pb.bam --sort --preset CCS --sample sample1 --rg '@RG\tID:movie1'
	"""
}

process ALIGN_ONT {
	tag "Create align for Oxord Nanopore Technology Callers"
	
	publishDir "${params.outdir}/aligment/${fastq.simpleName}", mode: "copy"
	
	input:
	path fastq
	
	output:
	path 'output_ont.bam', emit: bam
	path 'output_ont.bam.bai', emit: index
	val fastq.simpleName, emit: sample
	
	script:
	"""
	minimap2 -t ${params.threads} -ax map-ont --MD ${params.ref} $fastq > file.sam
	samtools view -@ ${params.threads} -S -b file.sam > file.bam
	samtools sort -@ ${params.threads} file.bam -o file.sorted.bam
	samtools view -@ ${params.threads} -b -F 4 file.sorted.bam > output_ont.bam
	samtools index -@ ${params.threads} output_ont.bam
	rm file.sam file.bam file.sorted.bam
	"""
}

process PBSV {
	tag "Calling PBSV"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "pbsv.vcf"

	script:
	"""
	pbsv discover $bam pbsv.svsig.gz
	tabix -c '#' -s 3 -b 4 -e 4 pbsv.svsig.gz
	pbsv call -j ${params.threads} ${params.ref} pbsv.svsig.gz pbsv.vcf
	"""
}

process Sniffles {
	tag "Calling Sniffles"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "sniffles.vcf"

	script:
	"""
	sniffles -t ${params.threads} -i $bam -v sniffles.vcf --minsvlen 50
	"""
}

process CuteSV {
	tag "Calling CuteSV"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
	
	output:
	path "cuteSV.vcf"

	script:
	"""
	cuteSV -t ${params.threads} -l 50 -s 5 $bam ${params.ref} cuteSV.vcf .	
	"""
}

process Svim {
	tag "Calling Svim"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "svim.vcf"

	script:
	"""
	svim alignment . $bam ${params.ref}
	bcftools view -i 'QUAL >= 10' variants.vcf
	mv variants.vcf svim.vcf
	"""
}

process Dysgu {
	tag "Calling Dysgu"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "dysgu.vcf"
	
	script:
	"""
	dysgu call --mode nanopore ${params.ref} temp $bam > dysgu.vcf	
	"""
}

process Nanovar {
	tag "Calling Nanovar"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "nanovar.vcf"

	script:
	"""
	nanovar -t ${params.threads} -x ont $bam -l 50 ${params.ref} . --mdb /tools/ncbi-blast-2.3.0+/bin/makeblastdb --wmk /tools//ncbi-blast-2.3.0+/bin/windowmasker --hsb /tools/queries/hs-blastn-src/v0.0.5/hs-blastn
	mv output_ont.nanovar.pass.vcf nanovar.vcf
	"""
}

process NanoSV {
	tag "Calling NanoSV"

	publishDir "${params.outdir}/vcfs/${sample}", mode: "copy"

	input:
	path bam
	val sample
	path bai
	
	output:
	path "nanoSV.vcf"

	script:
	"""
	NanoSV -t ${params.threads} -s /tools/samtools-1.12/samtools -c /tools/ConsensusSV-ONT-pipeline/config.ini -b /tools/ConsensusSV-ONT-pipeline/random_positions_chrxy.bed $bam -o nanoSV.vcf
	"""
}

