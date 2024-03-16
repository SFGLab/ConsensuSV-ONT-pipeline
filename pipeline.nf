#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and fastq files and can be provided as command line options
 */

params.input = "input.csv"	
params.outdir = "results"
params.ref = "/tools/GRCh38_full_analysis_set_plus_decoy_hla.fa"
params.threads = 50
params.mem = 500

workflow {
    	files = Channel.fromPath(params.input).splitCsv()

	ALIGN_PB(files)
	ALIGN_ONT(files)
	
	PBSV(ALIGN_PB.out.bam, ALIGN_PB.out.sample, ALIGN_PB.out.index)
	Sniffles(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	CuteSV(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Svim(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Dysgu(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	Nanovar(ALIGN_ONT.out.bam, ALIGN_ONT.out.sample, ALIGN_ONT.out.index)
	
	ConsensuSV_ONT(ALIGN_ONT.out.sample, PBSV.out.vcf_filtered, Sniffles.out.vcf_filtered, CuteSV.out.vcf_filtered, Svim.out.vcf_filtered, Dysgu.out.vcf_filtered, Nanovar.out.vcf_filtered, ALIGN_ONT.out.bam, ALIGN_ONT.out.index)

}

process ALIGN_PB {
	tag "Create align for Pac Bio Callers"

	publishDir "${params.outdir}/aligment/${fastq.simpleName}"
	
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
	
	publishDir "${params.outdir}/aligment/${fastq.simpleName}"
	
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

	publishDir "${params.outdir}/vcfs/${sample}"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "pbsv.vcf", emit: vcf
	path "pbsv_filtered.vcf", emit: vcf_filtered

	script:
	"""
	pbsv discover $bam pbsv.svsig.gz
	tabix -c '#' -s 3 -b 4 -e 4 pbsv.svsig.gz
	pbsv call -j ${params.threads} ${params.ref} pbsv.svsig.gz pbsv.vcf
	bcftools view -i 'INFO/IMPRECISE=0' pbsv.vcf > pbsv_filtered.vcf
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
	path "sniffles.vcf", emit: vcf
	path "sniffles_filtered.vcf", emit: vcf_filtered

	script:
	"""
	sniffles -t ${params.threads} -i $bam -v sniffles.vcf --minsvlen 50
	bcftools view -i 'INFO/IMPRECISE=0' sniffles.vcf > sniffles_filtered.vcf
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
	path "cuteSV.vcf", emit: vcf
	path "cuteSV_filtered.vcf", emit: vcf_filtered

	script:
	"""
	cuteSV -t ${params.threads} -l 50 -s 5 $bam ${params.ref} cuteSV.vcf .	
	bcftools view -i 'INFO/RE>10' cuteSV.vcf > cuteSV_filtered.vcf
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
	path "svim.vcf", emit: vcf
	path "svim_filtered.vcf", emit: vcf_filtered

	script:
	"""
	svim alignment . $bam ${params.ref}
	bcftools sort variants.vcf > svim.vcf
	bcftools view -i 'QUAL>=10' svim.vcf > svim_filtered.vcf
	"""
}

process Dysgu {
	tag "Calling Dysgu"

	publishDir "${params.outdir}/vcfs/${sample}"

	input:
	path bam
	val sample
	path bai
 
	output:
	path "dysgu.vcf", emit: vcf
	path "dysgu_filtered.vcf", emit: vcf_filtered
	
	script:
	"""
	dysgu call --mode nanopore ${params.ref} temp $bam > dysgu.vcf	
	bcftools view -i 'FORMAT/PROB>0.5' dysgu.vcf > dysgu_filtered.vcf
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
	path "nanovar.vcf", emit: vcf
	path "nanovar_filtered.vcf", emit: vcf_filtered

	script:
	"""
	nanovar -t ${params.threads} -x ont $bam -l 50 ${params.ref} . --mdb /tools/ncbi-blast-2.3.0+/bin/makeblastdb --wmk /tools//ncbi-blast-2.3.0+/bin/windowmasker --hsb /tools/queries/hs-blastn-src/v0.0.5/hs-blastn
	mv output_ont.nanovar.pass.vcf nanovar.vcf
	bcftools view -i 'INFO/NN>0.5 & INFO/SVLEN!="."' nanovar.vcf > nanovar_filtered.vcf
	"""
}

process ConsensuSV_ONT {
	tag "Truvari variants merging and collapsing"

	publishDir "${params.outdir}/vcfs/${sample}"

	input:
	val sample
	path pbsv_vcf
	path sniffles_vcf
	path cutesv_vcf
	path svim_vcf
	path dysgu_vcf
	path nanovar_vcf
	path bam
	path index
	
	output:
	path 'consensuSV-ONT_DEL.vcf'
	path 'consensuSV-ONT_INS.vcf'

	script:
	"""
        bgzip -c  $nanovar_vcf > ${nanovar_vcf}.gz
        tabix -p vcf ${nanovar_vcf}.gz

	bgzip -c  $dysgu_vcf > ${dysgu_vcf}.gz
	tabix -p vcf ${dysgu_vcf}.gz
	
	bgzip -c  $svim_vcf > ${svim_vcf}.gz
	tabix -p vcf ${svim_vcf}.gz

	bgzip -c  $cutesv_vcf > ${cutesv_vcf}.gz
	tabix -p vcf ${cutesv_vcf}.gz

	bgzip -c  $sniffles_vcf > ${sniffles_vcf}.gz
	tabix -p vcf ${sniffles_vcf}.gz
	
	bgzip -c  $pbsv_vcf > ${pbsv_vcf}.gz
	tabix -p vcf ${pbsv_vcf}.gz
	
	bcftools merge -m none ${nanovar_vcf}.gz  ${pbsv_vcf}.gz ${sniffles_vcf}.gz ${cutesv_vcf}.gz ${svim_vcf}.gz ${dysgu_vcf}.gz --force-samples > ${sample}_merged.vcf
	
	bgzip -c  ${sample}_merged.vcf > ${sample}_merged.vcf.gz
	tabix -p vcf ${sample}_merged.vcf.gz

	bcftools view -i 'INFO/SVTYPE="DEL"' ${sample}_merged.vcf > ${sample}_merged_DEL.vcf
	bcftools view -i 'INFO/SVTYPE="INS" | INFO/SVTYPE="DUP"' ${sample}_merged.vcf >${sample}_merged_INS.vcf
	
	bgzip -c  ${sample}_merged_DEL.vcf > ${sample}_merged_DEL.vcf.gz
	tabix -p vcf ${sample}_merged_DEL.vcf.gz
	bgzip -c  ${sample}_merged_INS.vcf > ${sample}_merged_INS.vcf.gz
	tabix -p vcf ${sample}_merged_INS.vcf.gz
	
	truvari collapse -p=0 -i ${sample}_merged_DEL.vcf.gz -o ${sample}_truvari_merged_DEL.vcf -c ${sample}_truvari_collapsed_DEL.vcf
	truvari collapse -p=0 -i ${sample}_merged_INS.vcf.gz -o ${sample}_truvari_merged_INS.vcf -c ${sample}_truvari_collapsed_INS.vcf
	python /tools/ConsensusSV-ONT-pipeline/consensusv_ont.py -i ${sample}_truvari_merged_DEL.vcf -b $bam -t DEL -s $sample -o consensuSV-ONT_DEL.vcf
	python /tools/ConsensusSV-ONT-pipeline/consensusv_ont.py -i ${sample}_truvari_merged_INS.vcf -b $bam -t INS -s $sample -o consensuSV-ONT_INS.vcf
	
	"""
}

