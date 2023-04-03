#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and read pairs and can be provided as command line options
 */

params.ref = "/tools/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.outdir = "results"
params.design = "/tools/files.csv"	
params.threads = 4
params.mem = 4

workflow {
    files = Channel.fromPath(params.design).splitCsv()
    BAMFILE(files)
    INDEX(BAMFILE.out.bam)
	Sniffles(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
	CuteSV(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
	Svim(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
	Dysgu(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
	//Nanovar(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
	NanoSV(BAMFILE.out.bam, BAMFILE.out.sample, INDEX.out)
}

process BAMFILE {
    tag "Create symlink to bam files"
 
    input:
    path bam
 
    output:
    path 'output.bam', emit: bam
    val bam.simpleName, emit:sample
 
    script:
    """
	ln -s $bam output.bam
    """
}
 

process INDEX {
    tag "Indexing files"
 
    input:
    path bam
 
    output:
    path "${bam}.bai"
 
    script:
    """
    samtools index $bam
    """
} 

process Sniffles {
    tag "Calling Sniffles"

    publishDir "${params.outdir}/vcfs/${sample}"

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

    publishDir "${params.outdir}/vcfs/${sample}"

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

    publishDir "${params.outdir}/vcfs/${sample}"

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
    tag "Calling Svim"

    publishDir "${params.outdir}/vcfs/${sample}"

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

    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "nanovar.vcf"

    script:
    """
	nanovar -t ${params.threads} -x ont $bam -l 50 ${params.ref} . --mdb /tools/ncbi-blast-2.3.0+/bin/makeblastdb --wmk /tools//ncbi-blast-2.3.0+/bin/windowmasker --hsb /tools/queries/hs-blastn-src/v0.0.5/hs-blastn
	"""
}

process NanoSV {
    tag "Calling NanoSV"

    publishDir "${params.outdir}/vcfs/${sample}"

    input:
    path bam
    val sample
    path bai
 
    output:
    path "nanoSV.vcf"

    script:
    """
	NanoSV -t ${params.threads} -s /tools/samtools-1.12/samtools -c /tools/config.ini -b /tools/random_positions_chrxy.bed $bam -o nanoSV.vcf
	"""
}

