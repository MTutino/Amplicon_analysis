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
   -r      Path to the directory with the reference databases, if not same as the script directory ***OPTIONAL***

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
export REF_DATA_PATH=

#Get the arguments
while getopts "hg:G:q:l:O:L:1:P:Si:o:m:t:r:3" OPTION
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
	 r)
	     if [ ! -d "$OPTARG" ] ; then
		 echo "Fatal: Non-existent directory '$OPTARG' supplied to -r option" >&2
		 exit 1
	     else
		 export REF_DATA_PATH=$(cd "$OPTARG" && pwd)
	     fi
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

# Get pipeline name in uppercase
PIPELINE_NAME=$(echo $PIPELINE | tr '[:lower:]' '[:upper:]')

#Check for required programs in current environment
#Report missing programs and exit if any are not installed/loaded
REQUIRED_PROGRAMS="cutadapt
sickle
fastqc
spades
bioawk
pandaseq
vsearch113
ChimeraSlayer.pl
print_qiime_config.py
fasta-splitter.pl
fasta_number.py
blastall
R"
# Require usearch executables for specific pipelines
case "$PIPELINE_NAME" in
    "UPARSE")
	# UPARSE pipeline needs usearch 8.0.1623
	REQUIRED_PROGRAMS="$REQUIRED_PROGRAMS
usearch8.0.1623_i86linux32"
	;;
    "QIIME")
	# UPARSE pipeline needs usearch 6.1.544
	REQUIRED_PROGRAMS="$REQUIRED_PROGRAMS
usearch6.1.544_i86linux32"
	;;
    *)
	;;
esac
MISSING_PROGRAMS=
for prog in $REQUIRED_PROGRAMS ; do
    echo -n "Checking for ${prog}..."
    if [ -z "$(which $prog)" ] ; then
	echo missing
	MISSING_PROGRAMS=yes
    else
	echo ok
    fi
done
if [ ! -z "$MISSING_PROGRAMS" ] ; then
    echo "One or more required programs are missing"
    exit 1
fi

#If $STEP3 is empty (the argument -3 has not been passed)  then the script will run the entire pipeline
if [[ -z $STEP3 ]]; then
		
		#Define the number of cores and the amount of RAM. If the ammount of ram is less than 22GB than do not use Silva.
		#Depending on the size of your file, less than 22GB might not be enought	
		if [ -z "$NSLOTS" ]; then
				export NSLOTS=$(nproc);
		fi;
		TOT_RAM=$(echo $(($(grep MemTotal /proc/meminfo | awk '{print $2}')/1024)));
		RAM_PER_CORE=$(($TOT_RAM/$NSLOTS));
		
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

		# Reference databases directory
		export REF_DATA_PATH="${REF_DATA_PATH:-$DIR}"
		echo "Expecting reference databases under $REF_DATA_PATH"
	
		
	
		# Define the path for every reference we are going to use. Not all of them are actually used.
		if [[ -z $SILVA ]]; then
				echo "Reference database: GreenGenes 13_8" >> $LOG;
				#Greengenes
				export REF="$REF_DATA_PATH/gg_13_8_otus/rep_set/97_otus.fasta";
				export TAX="$REF_DATA_PATH/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt";
				export TREE="$REF_DATA_PATH/gg_13_8_otus/trees/97_otus.tree";
				export ALIGNED="$REF_DATA_PATH/gg_13_8_otus/rep_set_aligned/97_otus.fasta";
				export CORE="$REF_DATA_PATH/gg_13_8_otus/rep_set_aligned/97_otus.fasta" ;
				export CHIM="$REF_DATA_PATH/RDPClassifier_16S_trainsetNo14_rawtrainingdata/trainset14_032015.fasta";
		else	
				echo "Reference database: Silva_119" >> $LOG;
				#Silva
				export REF="$REF_DATA_PATH/Silva/Silva119_release/rep_set/97/Silva_119_rep_set97.fna";
				export TAX="$REF_DATA_PATH/Silva/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt";
				export TREE="$REF_DATA_PATH/Silva/Silva119_release/97_FastTree_trees/Silva_119_rep_set97_aligned_16S_only_pfiltered.tre";
				export ALIGNED="$REF_DATA_PATH/Silva/Silva119_release_aligned_rep_files/97_16S_only/Silva_119_rep_set97_aligned_16S_only.fna";
				export CORE="$REF_DATA_PATH/Silva/Silva119_release/core_alignment/core_Silva119_alignment.fna";
				export CHIM="$REF_DATA_PATH/RDPClassifier_16S_trainsetNo14_rawtrainingdata/trainset14_032015.fasta";
		fi;
		# Check that the reference data actually exists
		echo "Checking for reference databases:"
		MISSING_DATABASES=
		for reference_db in $REF $TAX $TREE $ALIGNED $CORE $CHIM ; do
		    echo -n "${reference_db}..."
		    if [ ! -e $reference_db ] ; then
			echo missing
			MISSING_DATABASES=yes
		    else
			echo ok
		    fi
		done
		if [ ! -z "$MISSING_DATABASES" ] ; then
		    echo "One or more reference databases are missing" >&2
		    exit 1
		fi

################################################################## SECOND_STEP ###################################################################
		# Execute appropriate pipeline based on -P option
		case "$PIPELINE_NAME" in
		    "QIIME")
			source $DIR/QIIME.sh
			;;
		    "UPARSE")
			source $DIR/UPARSE.sh
			;;
		    "VSEARCH")
			source $DIR/VSEARCH.sh
			;;
		    *)
			# Unrecognised pipeline
			if [ -z $PIPELINE ]; then
			    #Exit if no pipeline specified
			    usage
			    exit 1
			else
			    #Invalid pipeline option
			    echo " -P $PIPELINE is not a valid option." >> $LOG;
			    usage;
			    exit 1;
			fi
			;;
		esac
################################################################ THIRD_STEP ####################################################################
		source $DIR/THIRD_STEP.sh

else
############################################################# ONLY_THIRD_STEP ####################################################################
		#If $STEP3 (-3) has been passed then only the last step will be run 
		source $DIR/STAT_ANALYSIS.sh STAT_ANALYSIS.sh
				
fi;
	
