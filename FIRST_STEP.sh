#!/bin/bash 
	
		echo "Quality control step started at $(date|awk '{print $4}')" >> $LOG;
		
				
		# Write parameters in a file
		echo -e "-q is:"$PHRED"\n""-l is:"$WINDOW"\n""-O is:"$OVERLAP"\n""-L is:"$LENGTH"\n""-g is:"$FORWARD"\n""-G is:"$REVERSE >> $LOG
		
		# If $NO_CUTADAPT is empty (forward and reverse passed) than Cutadapt will be used
		if [[ -z $NO_CUTADAPT ]]; then		
				echo "Cutadapt version 1.8.1. Started at $(date|awk '{print $4}')" >> $LOG;
			set +e
				for i in $(ls -d *|sed 's/\///g' 2> /dev/null);do 
						R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
						R2=$(ls $i/*_R2_*.fastq* 2> /dev/null);

						if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]];then 
								mkdir -p $i/Cutadapt;
								set -e;
								cutadapt -e 0.15 -g ^$FORWARD -G ^$REVERSE --untrimmed-output $i/Cutadapt/$(basename ${i})"_R1_untrimmed.fastq" --untrimmed-paired-output $i/Cutadapt/$(basename ${i})"_R2_untrimmed.fastq" -o $i/Cutadapt/$(basename ${i})"_R1_trimmed.fastq" -p $i/Cutadapt/$(basename ${i})"_R2_trimmed.fastq" $R1 $R2;
								set +e;
						fi;

				done;
			
		else
				echo "You did not pass any primers and Cutadapt has been disabled $(date|awk '{print $4}')" >> $LOG;
		fi;

		# Sickle
		echo "Sickle version 1.33. Started at $(date|awk '{print $4}')" >> $LOG;
		if [[ -z $NO_CUTADAPT ]]; then					
				for i in $(ls -d *|sed 's/\///g' 2> /dev/null);do 
						R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
						if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]]; then 
								mkdir -p $i/Cutadapt_Sickle/;
								mkdir -p $i/Cutadapt_Sickle/Q$PHRED/;
								set -e;
								sickle pe -f $i/Cutadapt/$(basename ${i})"_R1_trimmed.fastq" -r $i/Cutadapt/$(basename ${i})"_R2_trimmed.fastq" -o $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_R1_cutsick.fastq" -p $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_R2_cutsick.fastq" -s $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_singlet.fastq" -x -q $PHRED -l $WINDOW -t "sanger";
								set +e;
						fi;

				done;
			
		else
				set +e;
				for i in $(ls -d *|sed 's/\///g' 2> /dev/null); do 
						R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
						if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]]; then
								mkdir -p $i/Cutadapt_Sickle/;
								mkdir -p $i/Cutadapt_Sickle/Q$PHRED/; 
								set -e;
								sickle pe -f $i/*_R1_*.fastq* -r $i/*_R2_*.fastq* -o $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_R1_cutsick.fastq" -p $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_R2_cutsick.fastq" -s $i/Cutadapt_Sickle/Q$PHRED/$(basename ${i})"_singlet.fastq" -x -q $PHRED -l $WINDOW -t "sanger"; 
								set +e;
						fi;
				done;
		fi;

		#FastQC raw and trimmed data
		echo "FastQC version 0.11.3. Started at $(date|awk '{print $4}')" >> $LOG;

		for i in $(ls -d *|sed 's/\///g' 2> /dev/null); do 
				R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
				if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]]; then 
						mkdir -p $i/FastQC/;
						mkdir -p $i/FastQC/Raw;
						mkdir -p $i/FastQC/cutdapt_sickle/;
						mkdir -p $i/FastQC/cutdapt_sickle/Q$PHRED;
			
						#FastQC raw data
						R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
						R2=$(ls $i/*_R2_*.fastq* 2> /dev/null); 
						set -e;
						fastqc -o $i/FastQC/Raw $R1 $R2;
						set +e;	
						#FastQC trimmed data
						R3=$(ls $i/Cutadapt_Sickle/Q$PHRED/*_R1_*cutsick.fastq 2> /dev/null);
						R4=$(ls $i/Cutadapt_Sickle/Q$PHRED/*_R2_*cutsick.fastq 2> /dev/null);
						set -e;
						fastqc -o $i/FastQC/cutdapt_sickle/Q$PHRED $R3 $R4;
						set +e;
				fi;
		done;


		#SPAdes for error correction of illumina miseq reads
		echo "SPADes version 3.5.0. Started at $(date|awk '{print $4}')" >> $LOG;

		for i in $(ls -d *|sed 's/\///g' 2> /dev/null);do 
				R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
				if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]];then
						cd $i/Cutadapt_Sickle/Q$PHRED/; 
						set -e
						spades.py -1 *_R1_*cutsick.fastq -2 *_R2_*cutsick.fastq -o . --only-error-correction --careful --disable-gzip-output 2> /dev/null; 
						set +e;
						cd ../../../; 
				fi;
		done;
	
	
		# Pandaseq
		echo "Pandaseq version 2.8. Started at $(date|awk '{print $4}')" >> $LOG;

		for i in $(ls -d *|sed 's/\///g' 2> /dev/null);do 
				R1=$(ls $i/*_R1_*.fastq* 2> /dev/null); 
				if [[ $R1 == *_R1_*.fastq.gz ]] || [[ $R1 == *_R1_*.fastq ]];then  
					awk '/^@M01070/{$0=$0" 1:N:0:GGACTCCTGTAAGGAG"}1' $i/Cutadapt_Sickle/Q$PHRED/corrected/*_R1_*.cor.fastq > /tmp/$(basename ${i})"_forward.fastq";
					awk '/^@M01070/{$0=$0" 2:N:0:GGACTCCTGTAAGGAG"}1' $i/Cutadapt_Sickle/Q$PHRED/corrected/*_R2_*.cor.fastq > /tmp/$(basename ${i})"_reverse.fastq";
					set -e;
					pandaseq -f /tmp/$(basename ${i})"_forward.fastq" -r /tmp/$(basename ${i})"_reverse.fastq" -N -l $LENGTH -B -d bfsrk -o $OVERLAP > $i/$(basename ${i})".overlap.fasta" -g $i/pandaseq_log.txt; 
					set +e;
					rm /tmp/$(basename ${i})"_forward.fastq"; 
					rm /tmp/$(basename ${i})"_reverse.fastq"; 
				fi;
		done;
	
	
		# Stats from raw data to overlap sequences
		echo "Creating stats from raw data to overlap sequences $(date|awk '{print $4}')" >> $LOG;
		mkdir -p Stats ;
	
		echo -e "SAMPLE\tRaw_data\tPRIMER_CHECKED_PE_READS\tTRIMMED_PE_READS\tPE_READS_WITH_CHANGED_BASES\tCHANGED_BASES\tFAILED_BASES\tTOTAL_BASES\tFINAL_PE_READS\tOVERLAP_READS" >> Stats/Reads_count.txt ;

for i in $(ls -d */ 2> /dev/null); do 
				set -e;
				cd $i; 
				if [ -s $(basename ${i})".overlap.fasta" ]; then 
						(echo -ne "\n"$(basename ${i}); 
						if [ *_R1_*.fastq* == *_R1_*.fastq.gz ]; then
								echo -ne '\t'$(($(zcat *_R1_*.fastq.gz| wc -l)/4));
						else
								echo -ne '\t'$(($(cat *_R1_*.fastq| wc -l)/4));
						fi;
					
						if [[ -z $NO_CUTADAPT ]]; then
								echo -ne '\t'$(($(wc -l < Cutadapt/$(basename ${i})"_R1_trimmed.fastq")/4));
						fi;	
					
						echo -ne '\t'$(($(wc -l < Cutadapt_Sickle/Q$PHRED/*_R1_*cutsick.fastq)/4));
						echo -ne '\t'$(grep -i "Correction done" Cutadapt_Sickle/Q$PHRED/spades.log | grep -Po '(?<=in ).*(?= reads)'); 
						echo -ne '\t'$(grep -i "Correction done" Cutadapt_Sickle/Q$PHRED/spades.log | grep -Po '(?<=Changed ).*(?= bases)'); 
						echo -ne '\t'$(grep -i "Failed to correct" Cutadapt_Sickle/Q$PHRED/spades.log | grep -Po '(?<=correct ).*(?= bases)'); 
						echo -ne '\t'$(grep -i "Failed to correct" Cutadapt_Sickle/Q$PHRED/spades.log | grep -Po '(?<=out of ).*(?=\.)'); 
						echo -ne '\t'$(($(wc -l < Cutadapt_Sickle/Q$PHRED/corrected/*_R1_*.cor.fastq)/4));
						echo -e '\t'$(grep -c ">" *.overlap.fasta) 
						)  
				fi >> ../Stats/Reads_count.txt; 
				cd ..; 
				set +e;
		done;
	
		#Length distribution for raw and trimmed data
		rm -f QUALITY_CONTROL/Length_distribution_raw.txt;
		rm -f QUALITY_CONTROL/Length_distribution_trimmed.txt;
		rm -f Stats/Length_distribution_raw.txt;
		rm -f Stats/Length_distribution_trimmed.txt;
		for i in $(ls -d */ 2> /dev/null); do
				set -e;
				if [[ -s $i/$(basename ${i})".overlap.fasta" ]]; then
						if [ $i/*_R1_*.fastq* == $i/*_R1_*.fastq.gz ]; then
								zcat $i/*_R1_*.fastq.gz >> Stats/total_reads_R1.fastq;
								zcat $i/*_R2_*.fastq.gz >> Stats/total_reads_R2.fastq;
						else
								cat $i/*_R1_*.fastq >> Stats/total_reads_R1.fastq;
								cat $i/*_R2_*.fastq >> Stats/total_reads_R2.fastq;
						fi;
						cat $i/Cutadapt_Sickle/Q$PHRED/*_R1_*cutsick.fastq >> Stats/total_reads_R1_trimmed.fastq;
						cat $i/Cutadapt_Sickle/Q$PHRED/*_R2_*cutsick.fastq >> Stats/total_reads_R2_trimmed.fastq;
				fi;
				set +e;
		done;
		set -e;
		bioawk -cfastx '{i[length($seq)]++}END{for(j in i)print j","i[j]}' Stats/total_reads_R1.fastq | sort -nrk2 -t","  >> Stats/Length_R1.txt;
		bioawk -cfastx '{i[length($seq)]++}END{for(j in i)print j","i[j]}' Stats/total_reads_R2.fastq | sort -nrk2 -t","  >> Stats/Length_R2.txt;
		bioawk -cfastx '{i[length($seq)]++}END{for(j in i)print j","i[j]}' Stats/total_reads_R1_trimmed.fastq | sort -nrk2 -t","  >> Stats/Length_R1_trimmed.txt;
		bioawk -cfastx '{i[length($seq)]++}END{for(j in i)print j","i[j]}' Stats/total_reads_R2_trimmed.fastq | sort -nrk2 -t","  >> Stats/Length_R2_trimmed.txt;
	
		((echo -e  $(column -t Stats/Length_R1.txt); 
		  echo -e $(column -t Stats/Length_R2.txt); 
		)|
		python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))"| 
		awk -F "," 'BEGIN {print "\t\t""Raw_data""\n""\t""Forward""\t\t""Reverse""\n""Length""\t\t""Count"" ""Length""\t\t""Count"} {print $1"\t\t"$2"\t\t"$3"\t\t"$4"\t\t"$5}'|
		expand --tabs=14
		) >> Stats/Length_distribution_raw.txt ;
		
		((echo -e $(column -t Stats/Length_R1_trimmed.txt); 
		  echo -e  $(column -t Stats/Length_R2_trimmed.txt)
		)|
		python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))"| 
		awk -F "," 'BEGIN {print "\t\t""Trimmed_data""\n""\t""Forward""\t\t""Reverse""\n""Length""\t\t""Count"" ""Length""\t\t""Count"} {print $1"\t\t"$2"\t\t"$3"\t\t"$4"\t\t"$5}'|
		expand --tabs=14
		) >> Stats/Length_distribution_trimmed.txt ;
	
	
	
	
		rm -f Stats/total_reads_R1.fastq;
		rm -f Stats/total_reads_R2.fastq;
		rm -f Stats/total_reads_R1_trimmed.fastq;
		rm -f Stats/total_reads_R2_trimmed.fastq;
		rm -f Stats/Length_R1.txt;
		rm -f Stats/Length_R2.txt;
		rm -f Stats/Length_R1_trimmed.txt;
		rm -f Stats/Length_R2_trimmed.txt;
		
		#Lets reorganize the outputs of the first step
		mkdir -p QUALITY_CONTROL ;
		rsync -a Stats/ QUALITY_CONTROL ;
		rm -r Stats/ ;
		
		##R plot needs to be optimised 
		#Create the Length distribution files with the right format for R
		#cd QUALITY_CONTROL;
		#tail -n+3 Length_distribution_raw.txt > Length_distribution_for_R_plots_raw.txt;
		#tail -n+3 Length_distribution_trimmed.txt > Length_distribution_for_R_plots_trimmed.txt;
		#Create plots with R.2 different scripts for version 2 and version 3 chemistry
		#if [ $(awk NR==2'{ print $1 }' Length_distribution_for_R_plots_raw.txt) == 251 ];then 
		#		Rscript $DIR/R_plots_v2.r;
		#elif  [ $(awk NR==2'{ print $1 }' Length_distribution_for_R_plots_raw.txt) == 301 ];then
		#		Rscript $DIR/R_plots_v3.r;
		#fi;
		
		#rm Length_distribution_for_R_plots_raw.txt;
		#rm Length_distribution_for_R_plots_trimmed.txt;
		#cd ..;
		
				
		echo "Quality control step finished at $(date|awk '{print $4}')" >> $LOG;
