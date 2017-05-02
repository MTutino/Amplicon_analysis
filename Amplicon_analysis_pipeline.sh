#!/bin/bash

#*****************************************************************************************************
#
# Name:			Amplicon_analysis_pipeline.sh
# Author:		Mauro Tutino (mauro.tutino@manchester.ac.uk)
# Purpose:		Analysis of 16S rRNA data from Illumina Miseq (Casava >= 1.8) paired-end reads.
#
#*****************************************************************************************************

# Best options for amplicon sequencing illumina miseq analysis V3-V4: -q 20 -l 10 -o 10 -L 380

set -o errexit

#Function usage
usage()
{
cat << EOF
usage: $0 options

######################################################################################################################################################################################
#
#   This script is for the analysis of 16S rRNA data from Illumina Miseq (Casava >= 1.8) paired-end reads. 
#	It is divided in 3 steps that can be run separately:
#
#	1)  Quality control with Cutadapt to check PCR primers, Sickle for quality trimming, SPADes for errors correction and Pandaseq to merge paired-end reads
#	2)  Singletons and chimeras removal, OTU clustering
#	
#		QIIME pipeline: Open-reference method with usearch61 if possible (it depends on the size), otherwise closed-reference method with uclust
#	
#		UPARSE pipeline: Free version cut off at 4GB RAM
#	
#		Vsearch pipeline: A freely available programme almost identical to Uparse but not limited to 4GB RAM
#
#	
#	3)  QIIME to perform beta and alpha diversity analysis
#
#
######################################################################################################################################################################################

The script requires the files "Final_name.txt" and "Metatable.txt".
"Categories.txt" is optional.

For an example of these files look at the README.txt.

Best options for the analysis of V3-V4 hypervariable regions: -q 20 -l 10 -o 10 -L 380
OPTIONS:
   -h      Show this message
   -g      Forward primer (Cutadapt) ***REQUIRED IF USING CUTADAPT. IF IT IS NOT PASSED, CUTADAPT WILL BE DISABLED***
   -G      Reverse primer (Cutadapt) ***REQUIRED IF USING CUTADAPT.IF IT IS NOT PASSED, CUTADAPT WILL BE DISABLED***
   -q      Phred score threshold below which the reads will be trimmed [default 20] (Sickle) ***OPTIONAL*** 
   -l (Lowercase "L")	   Length of the sliding Window in bp [default 10] (Sickle) ***OPTIONAL***
   -O (Uppercase "O")     Minimum overlap in bp between forward and reverse reads [default 10] (Pandaseq) ***OPTIONAL***
   -L      Minimum length in bp for a sequence to be kept after overlapping [default 380] (Pandaseq) ***OPTIONAL***
   -1 (One)     Use this option "-1 suppress" to skip the QC step
   -P	   Use this option to decide which pipeline to use, UPARSE, Vsearch or QIIME. UPARSE="-P uparse". Vsearch="-P vsearch". QIIME="-P QIIME"  ***REQUIRED***
   -S	   The default reference database is GreenGenes. Use this option without any argument if you want to use Silva. To use Silva you need at least 22 Gb of RAM.

 *** To run only the third step ***
   -3	   Pass this flag without any argument to run only the third step
   -i	   BIOM file
   -o (Lowercase "O")   Path to the directory where you want to store the result
   -m	   Metatable file. The names in the metatable has to match with those used to create the BIOM file
   -t	   Tree file. It has to be same as used in the pipeline (cannot use SILVA for the pipeline and GreenGenes for the last step)
   
   ***** IF YOU WANT TO USE SILVA DO NOT RUN THE SCRIPT WITH LESS THAN 5 CORES. ******
  *** If you want to modify the matatable file use Metatable_log/Metatable_mod.txt because sample names could be different, e.g. "-" instead of "_" ***
   
EXAMPLE USAGE:

## The following command will use default options for QC,vsearch pipeline and greengenes database
qsub Amplicon_analysis_pipeline.sh -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -P vsearch 

## The following command will NOT use cutadapt (for example to analyse samples amplified with different primers together) and it will use the uparse pipeline and silva database
qsub Amplicon_analysis_pipeline.sh -P uparse -S

## The following command will suppress the first step and the analysis will start directly from the second step.
## Useful if you have already analysed this data, for example with Vsearch, and you want to try a different one or if you used GreenGenes and you want to try Silva
qsub Amplicon_analysis_pipeline.sh -1 suppress -P uparse

## The following command will run only the third step
## Useful if you want to produce new plots after filtering out some samples from the OTU table.
qsub Amplicon_analysis_pipeline.sh -3 -i Path/Vsearch_OTU_tables/otu_table_without_mocks.biom -o Path/RESULTS_2/ -m Path/Metatable_log/Metatable_mod.txt -t Path/Vsearch_OTU_tables/otus.tre

################################################################################################################################################################################
#	IF YOU DID ANY MISTAKE OR SOMETHING WENT WRONG DURING THE FIRST STEP AND WANT TO START THE ANALYSIS FROM THE BEGINNING I ADVICE YOU TO DELETE EVERYTHING BEFORE TO START.
#	OTHERWISE SOME PROGRAMMES IN THE FIRST STEP COULD GET STUCK INSTEAD OF OVERWRITE THE FILES.
#################################################################################################################################################################################


EOF
};


#function check dimension
check_dimension()
	{
		# Check for multiplexed_lineazired.fasta dimension. Usearch61 32-bit is castrated to 4GB of RAM. For a file bigger then 3GB is likely to exceed this limit and it is better to use a closed-reference method.
		# Get dimension of the file
		export dimension=$(du -h Multiplexed_files/multiplexed_linearized.fasta|awk -F "\t" '{print $1}');

		echo "multiplexed_linearized.fasta is" $dimension >> $LOG;

		# The dimension could be in megabyte...
		regexm="^[0-9]+M$";
		if [[ "$dimension" =~ $regexm ]];then 
				export megabyte=$dimension;
		fi;

		# ...or in gigabyte
		regexgdot="^[0-9]+(\.[0-9]+)G$";
		regexg="^[0-9]+G$";
		if [[ "$dimension" =~ $regexg ]]; then
				gigabyte=$dimension;	
				export gigabyte=$(echo $gigabyte|sed 's/[^0-9]//g');
		elif [[ "$dimension" =~ $regexgdot ]]; then
				export gigabyte=$(echo $gigabyte|awk -F "." '{print $1}');
				
		fi;
		
	};
export -f check_dimension;
	
#Autoincrement log file. Everytime the script is launched it produces a different log file (log, log_1, log_2 etc)
export LOG=Amplicon_analysis_pipeline.log;
if [ -f $LOG ] && [ ! -f $LOG"_1" ]; then 
		export LOG=$LOG"_1";
fi;
	
if [ -f $LOG"_1" ];then 
		LOG=$LOG"_*";
		K=$(echo $LOG|awk '{print NF}');
		((++K));
		export LOG="Amplicon_analysis_pipeline.log_"$K;
fi;

		
#Initialise the variables
export FORWARD=
export REVERSE=
export PHRED=
export WINDOW=
export OVERLAP=
export LENGTH=
export STEP1=
export PIPELINE=
export NO_CUTADAPT=
export SILVA=
export MEMORY=
export STEP3=
export BIOM=
export RESULTS_PATH=
export METATABLE=
export TREE=

#Get the arguments
while getopts “hg:G:q:l:O:L:1:P:Si:o:m:t:3” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         g)
             export FORWARD=$OPTARG
             ;;
         G)
             export REVERSE=$OPTARG
             ;;
         q)
             export PHRED=$OPTARG
             ;;
         l)
             export WINDOW=$OPTARG
             ;;
         O)
             export OVERLAP=$OPTARG
             ;;
         L)
             export LENGTH=$OPTARG
             ;;
		 1)
             export STEP1=$OPTARG
             ;;
		 P)
             export PIPELINE=$OPTARG
             ;;
		 S)	 
			 export SILVA=silva
			 ;;
		 3)
             export STEP3=ONLY_LAST_STEP
             ;;
		 i)
             export BIOM=$OPTARG
             ;;
         o)
             export RESULTS_PATH=$OPTARG
             ;;
         m)
             export METATABLE=$OPTARG
             ;;
         t)
             export TREE=$OPTARG
             ;;
         \?)
             echo "invalid option: $OPTARG" >> $LOG
			 usage
             exit 1
             ;;
		 :)
			 echo "Option $OPTARG requires an argument." >> $LOG
			 exit 1
			 ;;
     esac
done;

#If no arguments are passed, exit the script
if [ -z "$1" ]; then
		usage
		exit 1
fi;

#Path to Script's directory
export DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
		
#Read in the list, create a variable for every programme and exit if the programme is not installed/loaded
LIST=$(echo -e "cutadapt\nsickle\nfastqc\nspades\nbioawk\npandaseq\nusearch8.0.1623_i86linux32\nvsearch113\nChimeraSlayer.pl\nusearch6.1.544_i86linux32\nprint_qiime_config.py\nfasta-splitter.pl\nfasta_number.py\nblastall\nR");
i=1;
while read -r programme; do 
		export PROG_$i="$programme";
		((++i));
done <<< "$LIST" ;

#I tried to use a while and a for loops but they did not work with "which"
PATH_1=$(which $PROG_1);
if [ -z "$PATH_1" ];then 
		echo "I could not find $PROG_1";
		echo "Make sure you can run $PROG_1 inside the working directory";
		exit 1;
fi;
PATH_2=$(which $PROG_2);
if [ -z "$PATH_2" ];then 
		echo "I could not find $PROG_2";
		echo "Make sure you can run $PROG_2 inside the working directory";
		exit 1;
fi;
PATH_3=$(which $PROG_3);
if [ -z "$PATH_3" ];then 
		echo "I could not find $PROG_3";
		echo "Make sure you can run $PROG_3 inside the working directory";
		exit 1;
fi;
PATH_4=$(which $PROG_4);
if [ -z "$PATH_4" ];then 
		echo "I could not find $PROG_4";
		echo "Make sure you can run $PROG_4 inside the working directory";
		exit 1;
fi;
PATH_5=$(which $PROG_5);
if [ -z "$PATH_5" ];then 
		echo "I could not find $PROG_5";
		echo "Make sure you can run $PROG_5 inside the working directory";
		exit 1;
fi;
PATH_6=$(which $PROG_6);
if [ -z "$PATH_6" ];then 
		echo "I could not find $PROG_6";
		echo "Make sure you can run $PR_6 inside the working directory";
		exit 1;
fi;
PATH_7=$(which $PROG_7);
if [ -z "$PATH_7" ];then 
		echo "I could not find $PROG_7";
		echo "Make sure you can run $PROG_7 inside the working directory";
		exit 1;
fi;
PATH_8=$(which $PROG_8);
if [ -z "$PATH_8" ];then 
		echo "I could not find $PROG_8";
		echo "Make sure you can run $PROG_8 inside the working directory";
		exit 1;
fi;
PATH_9=$(which $PROG_9);
if [ -z "$PATH_9" ];then 
		echo "I could not find $PROG_9";
		echo "Make sure you can run $PROG_9 inside the working directory";
		exit 1;
fi;
PATH_10=$(which $PROG_10);
if [ -z "$PATH_10" ];then 
		echo "I could not find $PROG_10";
		echo "Make sure you can run $PROG_10 inside the working directory";
		exit 1;
fi;
PATH_11=$(which $PROG_11);
if [ -z "$PATH_11" ];then 
		echo "I could not find $PROG_11";
		echo "Make sure you can run $PROG_11 inside the working directory";
		exit 1;
fi;
PATH_12=$(which $PROG_12);
if [ -z "$PATH_12" ];then 
		echo "I could not find $PROG_12";
		echo "Make sure you can run $PROG_12 inside the working directory";
		exit 1;
fi;
PATH_13=$(which $PROG_13);
if [ -z "$PATH_13" ];then 
		echo "I could not find $PROG_13";
		echo "Make sure you can run $PROG_13 inside the working directory";
		exit 1;
fi;
PATH_14=$(which $PROG_14);
if [ -z "$PATH_14" ];then 
		echo "I could not find $PROG_14";
		echo "Make sure you can run $PROG_14 inside the working directory";
		exit 1;
fi;
PATH_15=$(which $PROG_15);
if [ -z "$PATH_15" ];then 
		echo "I could not find $PROG_15";
		echo "Make sure you can run $PROG_15 inside the working directory";
		exit 1;
fi;

#If $STEP3 is empty (the argument -3 has not been passed)  then the script will run the entire pipeline
if [[ -z $STEP3 ]]; then
		
		#Define the number of cores and the amount of RAM. If the ammount of ram is less than 22GB than do not use Silva.
		#Depending on the size of your file, less than 22GB might not be enought	
		if [ -z "$NSLOTS" ]; then
				export NSLOTS=$(nproc);
				TOT_RAM=$(echo $(($(grep MemTotal /proc/meminfo | awk '{print $2}')/1024)));
				RAM_PER_CORE=$(($TOT_RAM/$NSLOTS));
		else
				RAM_PER_CORE="4000";
		fi;
		
		if [ $(($RAM_PER_CORE*$NSLOTS)) -lt 22000 ] && [[ -n $SILVA ]]; then
				echo -e "You need at least 22 gigabytes to use Silva\n The script will use GreenGenes instead\n At the end of the analysis, launch the script again with more cores." >> $LOG
				export SILVA=
		fi;

		echo "Job started at $(date|awk '{print $4}')" >> $LOG;

		# Check for parameters, if empty use the default values
		if [[ -z $PHRED ]]; then
				export PHRED="20" ;
		fi;

		if [[ -z $WINDOW ]]; then
				export WINDOW="10" ;
		fi;

		if [[ -z $OVERLAP ]]; then
				export OVERLAP="10" ;
		fi;

		if [[ -z $LENGTH ]]; then
				export LENGTH="380" ;
		fi;

		if [[ -z $SILVA ]]; then
				export MEMORY="4000";
		else	
				#export MEMORY="18000";
				export MEMORY="22000";
		fi;

		# Check for primers, if empty do not use Cutadapt
		if [[ -z $FORWARD ]] || [[ -z $REVERSE ]]; then
				export NO_CUTADAPT=no_cutadapt
		fi;



		#Check that final_name exists
		NAMES="Final_name.txt"
		if [ -f $NAMES ]; then
				echo "Checking for $NAMES" >> $LOG;
				echo "$NAMES exists" >> $LOG;
				sed -i 's/\r//g' $NAMES;
		else
				echo -e "I could not find the file \"$NAMES\". It is misspelled or it does not exist. This file is required to start the analysis." >> $LOG;
				exit 1 ;
		fi;

		#Check that there is a unique name per sample and that the samples are spelt right
		#Get the number of unique final names
		FIN_NAME=$(awk '{print $2}' $NAMES|awk NF|sort|uniq|wc -l);
		#Get the number of initial names (fastq files)
		INIT_NAME=$(awk '{print $1}' $NAMES|awk NF|wc -l);
		#There should not be identical names in the first column because of forward and reverse reads.
		UNIQUE=$(awk '{print $1}' $NAMES|sort|uniq|awk NF|wc -l);
		#If there are no unique names in the file then exit 
		if [[ $INIT_NAME -ne $UNIQUE ]]; then 
				echo "I found no unique names in the first column of $NAMES" >> $LOG;
				echo "A unique name per sample is required" >> $LOG;
				exit 1;
		fi;
		#If the number of initial names is not divisible by the number of final names that means there are not unique names. There should be a Fw and Rv for every final name
		if (($INIT_NAME%$FIN_NAME));then
				echo "I found no unique names in the second column of $NAMES" >> $LOG;
				echo "A unique name per sample is required" >> $LOG;
				exit 1;
		fi;

		#Check for the metatable. Without this file the analysis cannot start.
		METATABLE_1="Metatable.txt" ;
		echo "Checking for $METATABLE_1" >> $LOG;
		if [ -f $METATABLE_1 ]; then	 	
				echo "$METATABLE_1 exists" >> $LOG;
				sed -i 's/\r//g' $METATABLE_1;
				validate_mapping_file.py -m $METATABLE_1 -B -p -o Metatable_log/;
				#Check if there are errors in the Metatable.log file. If so, exit. Do not look at the warnings
				ERROR=$(awk 'NR==3 {print $1}' Metatable_log/Metatable.log);
				if [[ $ERROR == Found ]]; then 
						echo -e "Some errors have been found in the Metatable.txt file. Correct them before to run again the script.\nYou can find information about the errors in Metatable_log/Metatable.log file." >> $LOG;
						exit 1;
				fi;
				
				rm -f Metatable_log/Metatable_mod.txt;
				grep \# $METATABLE_1 > Metatable_log/Metatable_mod.txt;
				awk '{gsub(/_/,"-",$1);print}' $METATABLE_1|sed s/"\t""\t"/" "/g|sed s/" "/"\t"/g|grep -v \#|sed s/"\t"/"\t\t"/2 >> Metatable_log/Metatable_mod.txt ;
				export METATABLE="Metatable_log/Metatable_mod.txt" ;
		else
				echo "I could not find the file $METATABLE_1. It is misspelled or it does not exist."  >> $LOG;
				exit 1 ;
		fi;

		#Check if names in the first column of metatable and second column of Final_name match
		awk '{print $2}' $NAMES > Metatable_temp; grep -v \# $METATABLE_1|awk '{print $1}' >> Metatable_temp ;
		MERGED=$(echo "$(sort Metatable_temp|uniq|awk NF|wc -l)");
		META=$(grep -v \# $METATABLE_1|awk '{print $1}'|awk NF|wc -l);
		if [[ $MERGED -ne $META ]]; then 
				echo "Names in the first column of $METATABLE_1 and in the second column of $NAMES do not match" >> $LOG;
				exit 1;
		fi;
		rm -f Metatable_temp ;

		#Command to create folders from a list of file's names and folder's names. Folder's name will be the final name.
		set +e;	
		for i in $(ls); do 
				set -e;
				while read -r file_name final_name; do 
						if [[ "$i" == "$file_name" ]]; then 
								mkdir -p $final_name; 
								mv $i $final_name;
						fi;
				done < $NAMES;
				set +e;
		done;
		set -e;
################################################################# FIRST_STEP ###################################################################

		# Check for skipping quality step. If $STEP1 is empty the script will perform the quality control step
		if [[ -z $STEP1 ]]; then
				source $DIR/FIRST_STEP.sh
		else
				echo -e "You passed the option \"-1 suppress\" to disable the quality control step." >> $LOG;
		fi;

####################################################### JUNCTION BETWEEN FIRST AND SECOND STEP #################################################
		# Check if multiplexed and multiplexed_linearized.fasta already exist and if their size is not zero
		FILENAME1="Multiplexed_files/multiplexed.fasta";
		FILENAME2="Multiplexed_files/multiplexed_linearized.fasta";

		if [ -s $FILENAME1 ]; then

				echo "Multiplexed.fasta already exists. I will go straight to step 2." >> $LOG;
		else
				#The variable will contain the first and second column of metatable file (folder name and barcode)
				BARCODE=$(awk '{if (NR==1) next;print $1"\t"$2}' $METATABLE_1);
					
		
				#New and easier method to multiplex fasta files
				mkdir -p Multiplexed_files;
set +e
				for i in $(ls -d */|sed 's|[/]||g'); do 
set -e
						while read -r foldername barcode ; do
								if [[ "$i" == "$foldername" ]]; then 
										foldernamemod=$(echo $foldername|sed 's/_/-/g');
										awk -v k=$foldernamemod"_"$barcode '/^>/{gsub(">","",$0);split(k,a,"_");$0=">"a[1]"_"(++i)" "$0" orig_bc="a[2]" new_bc="a[2]" bc_diffs=0"}1' < $i/*overlap.fasta
								fi >> Multiplexed_files/multiplexed.fasta;
								oldnamemod=
						done <<< "$BARCODE";
				done;

				# Linearize multiplexed.fasta
				awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {print}' Multiplexed_files/multiplexed.fasta > Multiplexed_files/multiplexed_linearized.fasta;

		fi;
	
		
	
		# Define the path for every reference we are going to use. Not all of them are actually used.
		if [[ -z $SILVA ]]; then
				echo "Reference database: GreenGenes 13_8" >> $LOG;
				#Greengenes
				export REF="$DIR/gg_13_8_otus/rep_set/97_otus.fasta";
				export TAX="$DIR/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt";
				export TREE="$DIR/gg_13_8_otus/trees/97_otus.tree";
				export ALIGNED="$DIR/gg_13_8_otus/rep_set_aligned/97_otus.fasta";
				export CORE="$DIR/gg_13_8_otus/rep_set_aligned/97_otus.fasta" ;
				export CHIM="$DIR/RDPClassifier_16S_trainsetNo14_rawtrainingdata/trainset14_032015.fasta";
		else	
				echo "Reference database: Silva_119" >> $LOG;
				#Silva
				export REF="$DIR/Silva/Silva119_release/rep_set/97/Silva_119_rep_set97.fna";
				export TAX="$DIR/Silva/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt";
				export TREE="$DIR/Silva/Silva119_release/97_FastTree_trees/Silva_119_rep_set97_aligned_16S_only_pfiltered.tre";
				export ALIGNED="$DIR/Silva/Silva119_release_aligned_rep_files/97_16S_only/Silva_119_rep_set97_aligned_16S_only.fna";
				export CORE="$DIR/Silva/Silva119_release/core_alignment/core_Silva119_alignment.fna";
				export CHIM="$DIR/RDPClassifier_16S_trainsetNo14_rawtrainingdata/trainset14_032015.fasta";
		fi;
################################################################## SECOND_STEP ###################################################################
		#Exit if no pipeline specified
		if [ -z $PIPELINE ]; then
				usage;
				exit 1;
		fi;
		# If -P is passed as "qiime" the script will use QIIME
		if [[ "$PIPELINE" == "qiime" ]] || [[ "$PIPELINE" == "QIIME" ]] || [[ "$PIPELINE" == "Qiime" ]]; then		
				source $DIR/QIIME.sh
		
		# If -P is passed as "uparse" the script will use UPARSE (Usearch 8.0)
		elif [[ "$PIPELINE" == "UPARSE" ]] || [[ "$PIPELINE" == "uparse" ]] || [[ "$PIPELINE" == "Uparse" ]]; then
				source $DIR/UPARSE.sh
		
		# If -P is passed as "vsearch" the script will use vsearch, a freely available programme (64-bit) almost identical to UPARSE		
		elif  [[ "$PIPELINE" == "VSEARCH" ]] || [[ "$PIPELINE" == "vsearch" ]] || [[ "$PIPELINE" == "Vsearch" ]]; then	
				source $DIR/VSEARCH.sh
		
		else
				echo " -P $PIPELINE is not a valid option." >> $LOG;
				usage;
				exit 1;		
		fi;
################################################################ THIRD_STEP ####################################################################
		source $DIR/THIRD_STEP.sh

else
############################################################# ONLY_THIRD_STEP ####################################################################
		#If $STEP3 (-3) has been passed then only the last step will be run 
		source $DIR/STAT_ANALYSIS.sh STAT_ANALYSIS.sh
				
fi;
	
