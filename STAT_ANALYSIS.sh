#!/bin/bash

set -o errexit

#Function usage
usage()
{
cat << EOF
usage: Only the third step

After running Amplicon_analysis_pipeline.sh you may want to delete some samples from the analysis.
After doing that you can use this script to repeat just the last step and get new plots.
To use this script you need to pass the following arguments:

-i Path to the biom file

-o Path to a directory to store the results

-m Path to the metatable file

-t Path to the tree file

The tree file has to be the same you used for the analysis.
If you want to modify the matatable file use Metatable_log/Metatable_mod.txt because sample names could be different, e.g. "-" instead of "_"

Example:

## The following command will run only the third step
## Useful if you want to produce new plots after filtering out some samples from the OTU table.

Path_to_the_script/Amplicon_analysis_pipeline.sh -3 -i Path/Vsearch_OTU_tables/otu_table_without_mocks.biom -o Path/RESULTS_2/ -m Path/Metatable_log/Metatable_mod.txt -t Path/Vsearch_OTU_tables/otus.tre

EOF
};
		
		#Autoincrement log file. Everytime the script is launched it produces a new log file (log, log_1, log_2 etc)
		export LOG=Only_third_step.log;
		if [ -f $LOG ] && [ ! -f $LOG"_1" ]; then 
				export LOG=$LOG"_1";
		fi;
	
		if [ -f $LOG"_1" ];then 
				LOG=$LOG"_*";
				K=$(echo $LOG|awk '{print NF}');
				((++K));
				export LOG="Only_third_step.log_"$K;
		fi;
		
		echo "You passed the flag "-3" to run only the third step";
		#Chech if missing arguments. All of them are required
		if [[ -z $BIOM ]] || [[ -z $RESULTS_PATH ]] || [[ -z $METATABLE ]] || [[ -z $TREE ]]; then
				echo "Found missing arguments" >> $LOG;
				usage;
				exit 1;
		fi;
		
		echo "Check if files exist:"
		if [ -f $METATABLE ]; then	 	
				echo "$METATABLE exists" >> $LOG;
		else
				echo "I could not find the file $METATABLE. It is misspelled or it does not exist."  >> $LOG;
				exit 1 ;
		fi;
		if [ -f $BIOM ]; then	 	
				echo "$BIOM exists" >> $LOG;
		else
				echo "I could not find the file $BIOM. It is misspelled or it does not exist."  >> $LOG;
				exit 1 ;
		fi;

		if [ -f $TREE ]; then	 	
				echo "$TREE exists" >> $LOG;
		else
				echo "I could not find the file $TREE. It is misspelled or it does not exist."  >> $LOG;
				exit 1 ;
		fi;

		#Create the results directory if it does not exist
		mkdir -p $RESULTS_PATH;

		
		#Start the analysis
		source $DIR/THIRD_STEP.sh;
