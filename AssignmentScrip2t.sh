#!/bin/bash
		#PBS -l nodes=1:ppn=4:centos6,cput=24:00:00,walltime=48:00:00
		#PBS -N AssignmentScript
		#PBS -d /export/home/biostuds/2505621h/Assignment/Data
		#PBS -m abe
		#PBS -M 2505621h@student.gla.ac.uk
		#PBS -q bioinf-stud

	#DATA PATHS
		fqpath='/export/projects/polyomics/buzz/biostuds'
		adapter='/export/home/biostuds/2505621h/Assignment/Data/illumina_adapter.fa'
		hs2index='/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2'
		annot_gtf='/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf'


	#SUB DIRS
		qc_output_dir="/export/home/biostuds/2505621h/Assignment/FastQC_Results"
		mkdir -p $qc_output_dir
		hisat_dir='/export/home/biostuds/2505621h/Assignment/Hisat_Results'
		mkdir -p $hisat_dir
		stringtie_dir1='/export/home/biostuds/2505621h/Assignment/Stringtie_Results1'
		mkdir -p $stringtie_dir1


		for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12
		do
			raw_file="${sample}.c2.fq"
			link_name="${sample}.c2.fq"

			ln -s ${fqpath}/${raw_file} $link_name

		done


		for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12
		do
			fq_file="${sample}.c2.fq"

			/export/projects/polyomics/App/FastQC/fastqc -o $qc_output_dir $fq_file

		done


		for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12
		do
			input_file="${sample}.c2.fq"
			trim1="${sample}_tr1.fq"
			trim2="${sample}_tr2.fq"
			samfile="${sample}.sam"
			bamfile="${sample}.bam"
			bamfilesorted="${sample}.sort.bam"
			stringtiefile="${sample}.gtf"

			scythe -a illumina_adapter.fa $input_file -o $trim1 -q sanger

			sickle se -f $trim1 -o $trim2 -t sanger -q 10 -l 49

			hisat2 -p 4 --phred33 --rna-strandness R -x $hs2index -U $trim2 -S $hisat_dir/$samfile

			samtools view -bS -o $hisat_dir/$bamfile $hisat_dir/$samfile

			samtools sort -f $hisat_dir/$bamfile $hisat_dir/$bamfilesorted 

			rm $hisat_dir/$bamfile $hisat_dir/$samfile

			rm ${trim1} ${trim2}

			sampledir="${sample}_str_result1"
			strdir="${stringtie_dir1}/${sampledir}"

			mkdir -p $strdir

			stringtie $hisat_dir/$bamfilesorted --rf -p 4 -e -G $annot_gtf -B -o $strdir/$stringtiefile

		echo $sample $strdir/$stringtiefile >> samplelist.txt
		
		done
			

		python2.7 /export/projects/polyomics/App/prepDE.py	-i /export/home/biostuds/2505621h/Assignment/Data/samplelist.txt	









