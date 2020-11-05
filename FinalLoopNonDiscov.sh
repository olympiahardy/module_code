#!/bin/bash
		#PBS -l nodes=1:ppn=4:centos6,cput=24:00:00,walltime=48:00:00
		#PBS -N FinalScript
		#PBS -d /export/home/biostuds/2505621h/RNAseq/tx2dm/data
		#PBS -m abe
		#PBS -M 2505621h@student.gla.ac.uk
		#PBS -q bioinf-stud

		#RESOURCE FILES
		adapter='/export/home/biostuds/2505621h/RNAseq/tx2dm/data/illumina_adapter.fa'
		hs2index='/export/projects/polyomics/Genome/Dmelanogaster/dm3/Hisat2Index/dm3'
		gtf='/export/projects/polyomics/Genome/Dmelanogaster/dm3/dmel5.57/dmel-all-no-analysis-r5.57clean2.gtf'
		data='/export/home/biostuds/2505621h/RNAseq/tx2dm/data'

		# SUBDIRS
		hisat_dir='/export/home/biostuds/2505621h/RNAseq/tx2dm/hisat_results'
		stringtie_dir='/export/home/biostuds/2505621h/RNAseq/tx2dm/stringtie_results/' # ignore
		stringtie_dir2='/export/home/biostuds/2505621h/RNAseq/tx2dm/stringtie_results2/' # path to my directory for the stringtie results in the first loop
		stringtie_dir3='/export/home/biostuds/2505621h/RNAseq/tx2dm/stringtie_results3/' # path to my directory for the stringtie results in the second loop

		#RUNNING LOOP

		gtflist="" # empty initialisation of gtflist
		mergestr="stringtie_merged.gtf" # variable made for what I want my merged .gtf file to be called as a result of stringtie merge

		for sample in tb1 tb2 tb3 wf1 wf2 wf4
		do
			fastq="${sample}.fq"
			trim1="${sample}.t1.fq"
			trim2="${sample}.t2.fq"
			samfile="${sample}.sam"
			bamfile="${sample}.bam"
			bamfilesorted="${sample}.sort.bam"
			stringtiefile="${sample}.gtf"
			sampledir="${sample}_str_result2" # name of my sample specific directories I created on the command line (example script from the lab makes the directories in the loop but I found it easier to not do that!)

			scythe -a illumina_adapter.fa $fastq -o $trim1 -q illumina

			sickle se -f $trim1 -o $trim2 -t illumina -q 10 -l 51

			hisat2 -p 4 --phred64 -x $hs2index -U $trim2 -S $hisat_dir/$samfile

			samtools view -bS -o $hisat_dir/$bamfile $hisat_dir/$samfile

			samtools sort -m 1000000000 -f $hisat_dir/$bamfile $hisat_dir/$bamfilesorted 

		#	rm $hisat_dir/$bamfile $hisat_dir/$samfile
		#	rm ${trim1} ${trim2}

			strdir="${stringtie_dir2}${sampledir}" # complete path to my sample-specific directories that can be found within my stringtie_results2 directory


			stringtie $hisat_dir/$bamfilesorted -p 4 -G $gtf -o $strdir/$stringtiefile # output here is so each .gtf file is placed in the correct sample-specific directory
																						# NOTE: I got confused about the arguments in this earlier but it is correct because this is non-discovery script if you see lab instructions we had to remove -B and -e arguments

			gtflist="${gtflist} ${strdir}/${stringtiefile}" # adding the path of the stringtie .gtf file for each sample that is found in each sample-specific directory (see strdir variable comment!)

		done


		#MERGING sample-specific transcriptomes
		stringtie --merge -p 4 -o $stringtie_dir2$mergestr -G $gtf $gtflist # -o output here is simply path to my directory I want it to be put in + what I want the file to be named as
																			# -G will be the path to the reference .gtf file + the variable gtflist that has the path of all the .gtf files for the samples made in the loop

		for sample in tb1 tb2 tb3 wf1 wf2 wf4 # NOTE: For this loop I made another stringtie results directory (stringtie_dir3 variable) with additional sample-specific directories within it to keep it separate from the stringtie results for the first loop
		do
			sampledir2="${sample}_str_result3" # Declare necessary variables you need to use inside this second loop
			strdir2="${stringtie_dir3}${sampledir2}"
			bamfilesorted="${sample}.sort.bam"
			stringtiefile2="${sample}.gtf"

		 stringtie $hisat_dir/$bamfilesorted -p 4 -e -G $stringtie_dir2$mergestr -B -o $strdir2/$stringtiefile2 # input is sorted BAM files
		 																										# -G is the path to find my merged .gtf file
		 																										# output is path to sample specific directories like the stringtie output in the first loop 

		done



