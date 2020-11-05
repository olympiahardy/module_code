#!/bin/bash
		#PBS -l nodes=1:ppn=2:bio7,cput=12:00:00,walltime=12:00:00
		#PBS -N PP_Script
		#PBS -d /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq
		#PBS -m abe
		#PBS -M 2505621h@student.gla.ac.uk
		#PBS -q bioinf-stud

	#SUB DIRS
#		qc_output_dir="/export/home/biostuds/2505621h/pathogenPolyomicsData/FastQC_Results"
#		mkdir -p $qc_output_dir

		# bowtie_dir='/export/home/biostuds/2505621h/pathogenPolyomicsData/Bowtie2_Results'
		# mkdir -p $bowtie_dir


	#	for sample in LmexWT LmexAmpB
	#	do
	#		fq1_file="${sample}_1.fastq.gz"
	#		fq2_file="${sample}_2.fastq.gz"


	#		/export/projects/polyomics/App/FastQC/fastqc -o $qc_output_dir $fq1_file
	#		/export/projects/polyomics/App/FastQC/fastqc -o $qc_output_dir $fq2_file

	#	done


		 # for sample in LmexWT LmexAmpB
		 # do
		 # 	input_file1="${sample}_1.fastq.gz"
			# input_file2="${sample}_2.fastq.gz"

			# trim_galore --phred64 --paired --illumina $input_file1 $input_file2
		
		 # done

		# 	bowtie2-build -f /export/home/biostuds/2505621h/pathogenPolyomicsData/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa /export/home/biostuds/2505621h/pathogenPolyomicsData/Bowtie2_Results/Lmexicana

		#	bowtie2 -p 2 --very-sensitive --phred64 -x /export/home/biostuds/2505621h/pathogenPolyomicsData/Bowtie2_Results/Lmexicana -1 LmexWT_1_val_1.fq.gz -2 LmexWT_2_val_2.fq.gz 2> WT.log | samtools view -buS | samtools sort -T WT > WT.bam

		#	bowtie2 -p 2 --very-sensitive --phred64 -x /export/home/biostuds/2505621h/pathogenPolyomicsData/Bowtie2_Results/Lmexicana -1 LmexAmpB_1_val_1.fq.gz -2 LmexAmpB_2_val_2.fq.gz 2> AmpB.log | samtools view -buS | samtools sort -T AmpB > AmpB.bam

		# for sample in WT AmpB
		# do
		# 	bamfile="${sample}.bam"

		# 	samtools index -b $bamfile 

		# done

		# samtools faidx /export/home/biostuds/2505621h/pathogenPolyomicsData/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa

		# bamaddrg -b WT.bam -s WT -b AmpB.bam -s AmpB | freebayes -f /export/home/biostuds/2505621h/pathogenPolyomicsData/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa -p 2 --stdin > snp2.variant


		# bamCoverage -b WT.bam -of bedgraph -o WT_cov.bw
		# bamCoverage -b AmpB.bam -of bedgraph -o AmpB_cov.bw

		# vcffilter -f "QUAL > 20" -s /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/snp2.variant > snp_filtered2.variant

		# picard CreateSequenceDictionary R=/export/home/biostuds/2505621h/pathogenPolyomicsData/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa
		# java -jar $VTC SO -i data=/export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/snp_filtered.vcf -c exact -s AmpB_unique=c[data[AmpB]:data[WT]] -R /export/home/biostuds/2505621h/pathogenPolyomicsData/Reference/TriTrypDB-25_LmexicanaMHOMGT2001U1103.fa -o /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/AmpB_unique_snps.vcf
		# 


		# snpEff build -c /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/SnpEff.config -gff3 -v Lmex

		snpEff -Xmx4g -no-intron -no SPLICE_SITE_REGION -no SPLICE_SITE_DONOR -no SPLICE_SITE_ACCEPTOR -c /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/SnpEff.config Lmex /export/home/biostuds/2505621h/pathogenPolyomicsData/DNAseq/AmpB_unique_snps.vcf > AmpB_annotated_variants_final.vcf

		cat AmpB_annotated_variants_final.vcf | vcfEffOnePerLine.pl | SnpSift extractFields - CHROM POS REF ALT "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "GEN[*].GT" | grep 'HIGH\|MODERATE' > non-synonymous.txt
































