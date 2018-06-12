#!/bin/bash

set -o errexit

echo "Analysis step started at $(date|awk '{print $4}')" >> $LOG;
	
if [ -f $RESULTS_PATH/table_summary.txt ]; then
		rm $RESULTS_PATH/table_summary.txt
fi;
if [ -f $RESULTS_PATH/OTUs_count.txt ]; then
		rm $RESULTS_PATH/OTUs_count.txt
fi;
	
# Print a little summary of the otu_table. Number of sequences per sample and number of observations per sample.
biom summarize-table --qualitative -i $BIOM -o $RESULTS_PATH/OTUs_count.txt;
biom summarize-table -i $BIOM -o $RESULTS_PATH/table_summary.txt; 

#Check if there are samples with less than 5000 sequences. These will be excluded because it would be a problem for subsequent analyses
#In case samples with <5K seqs are found a new OTU table will be made
#Get a list of samples with more and less than 5000 sequences
LESSFIVEK=$(tail -n+17 $RESULTS_PATH/table_summary.txt |awk -F ":" '{if ( $2 <= 5000) print $1"\t"$2 }'|sed 's/ //g');
MOREFIVEK=$(tail -n+17 $RESULTS_PATH/table_summary.txt |awk -F ":" '{if ( $2 >= 5000) print $1 }'|sed 's/ //g');

#If there are samples with less than 5000 sequences, filter them out from the OTU table for subsequent analyses
if [[ -n $LESSFIVEK ]];then 
		
		echo -e "The following samples have less than 5000 sequences and have being excluded from the analysis:\nSample_Name\tNumber_of_Sequences\n$LESSFIVEK" >> $LOG; 
		echo -e "All the analises will be performed using the new filtered OTU table $(echo $BIOM|sed 's/otu.*biom//g')/otu_table_more_than_5K.biom" >> $LOG;
		#Make a file with a list of samples to keep (>5000 sequences)
		echo "$MOREFIVEK" > $RESULTS_PATH/Samples_with_more_than_5K_sequences;
		filter_samples_from_otu_table.py -i $BIOM -o $(echo $BIOM|sed 's/otu.*biom//g')/otu_table_more_than_5K.biom --sample_id_fp $RESULTS_PATH/Samples_with_more_than_5K_sequences;
		
		#Changing the path to the OTU table with more than 5K
		BIOM=$(echo $BIOM|sed 's/otu.*biom//g')/otu_table_more_than_5K.biom;
		
		#Create a new summary for the selected samples. Delete the old ones first
		if [ -f $RESULTS_PATH/table_summary.txt ]; then
				rm $RESULTS_PATH/table_summary.txt;
		fi;
		if [ -f $RESULTS_PATH/OTUs_count.txt ]; then
				rm $RESULTS_PATH/OTUs_count.txt;
		fi;
	
		# Print a little summary of the otu_table. Number of sequences per sample and number of observations per sample.
		biom summarize-table --qualitative -i $BIOM -o $RESULTS_PATH/OTUs_count.txt;
		biom summarize-table -i $BIOM -o $RESULTS_PATH/table_summary.txt; 
fi;
	
#Extract the minimum depth to use in beta-diversity
MINDEPTH=$(grep Min $RESULTS_PATH/table_summary.txt |awk -F " " '{print $2}'|awk -F "." '{print $1}');
	
#Beta-diversity analysis
echo "beta_diversity_through_plots.py $(date|awk '{print $4}')" >> $LOG;
beta_diversity_through_plots.py -o $RESULTS_PATH/beta_div_even/ -i $BIOM -m $METATABLE -t $TREE -e $MINDEPTH --suppress_emperor_plots -f ;
		
#Make 2D plots of the rarefied data from beta_diversity_through_plots.py
echo "make_2d_plots.py $(date|awk '{print $4}')" >> $LOG;
make_2d_plots.py -i $RESULTS_PATH/beta_div_even/weighted_unifrac_pc.txt -m $METATABLE -k black -o $RESULTS_PATH/beta_div_even/weighted_2d_plot --scree ;
make_2d_plots.py -i $RESULTS_PATH/beta_div_even/unweighted_unifrac_pc.txt -m $METATABLE -k black -o $RESULTS_PATH/beta_div_even/unweighted_2d_plot --scree ;
	
#Calculate alpha diversity
echo "Calculate alpha diversity $(date|awk '{print $4}')" >> $LOG
mkdir -p $RESULTS_PATH/Alpha_diversity/ ;

#CHAO1 and ACE are metrics based on singletons counts. They cannot be used because singletons are already filtered out from the analysis and the result would not be correct
ALPHA_METRICS="observed_species,shannon,simpson,fisher_alpha";
rm -rf $RESULTS_PATH/Rarefied_otu_tables/;
multiple_rarefactions.py -i $BIOM -m 200 -x $MINDEPTH -s 2000 -n 10 -o $RESULTS_PATH/Rarefied_otu_tables/ ;
rm -rf $RESULTS_PATH/Alpha_diversity/Alpha_rarefactions/;
mkdir -p $RESULTS_PATH/Alpha_diversity/Alpha_rarefactions/ ;
alpha_diversity.py -i $RESULTS_PATH/Rarefied_otu_tables/ -m $ALPHA_METRICS -o $RESULTS_PATH/Alpha_diversity/Alpha_rarefactions/ -t $TREE ;
collate_alpha.py -i $RESULTS_PATH/Alpha_diversity/Alpha_rarefactions/ -o $RESULTS_PATH/Alpha_diversity/collated_alpha/ ;
make_rarefaction_plots.py -i $RESULTS_PATH/Alpha_diversity/collated_alpha/ -m $METATABLE -o $RESULTS_PATH/Alpha_diversity/rarefaction_curves/;

#Compare Alpha diversity
CATEGORIES="Categories.txt";
if [[ -f $CATEGORIES ]]; then
		echo -e "Categories.txt exists.\nCompare alpha diversity $(date|awk '{print $4}')" >> $LOG
		sed -i 's/\r//g' $CATEGORIES;
		mkdir -p $RESULTS_PATH/Alpha_diversity/Alpha_diversity_boxplot/ ;
		CATEGORIES=$(awk '{print $0}' Categories.txt);CATEGORIES=$(echo $CATEGORIES|sed 's/ /,/g');
		METRICS_COLUMN=$(echo "$ALPHA_METRICS"|sed 's/,/\n/g');
		
		#for every metric run compare_alpha_diversity.py
		while IFS= read -r metric; do
				compare_alpha_diversity.py -i $RESULTS_PATH/Alpha_diversity/collated_alpha/"$metric".txt -m $METATABLE -c $CATEGORIES -o $RESULTS_PATH/Alpha_diversity/Alpha_diversity_boxplot/Categories_"$metric";
		done <<<"$METRICS_COLUMN" ;
fi;
	
#Taxonomic summaries of samples
echo "summarize_taxa_through_plots.py $(date|awk '{print $4}')" >> $LOG;
summarize_taxa.py -i $BIOM -o $RESULTS_PATH/taxa -L 2,3,4,5,6,7 -a;
plot_taxa_summary.py -i $RESULTS_PATH/taxa/*L2.txt,$RESULTS_PATH/taxa/*L6.txt,$RESULTS_PATH/taxa/*L7.txt -l Phylum,Genus,Species -c pie,bar -o $RESULTS_PATH/phylum_genus_charts/ ;
		

#Create the heatmap in html format		
make_otu_heatmap_html.py -i $BIOM -o $RESULTS_PATH/Heatmap -t $TREE  -m $METATABLE;	
	
echo "Third step finished at $(date|awk '{print $4}')" >> $LOG;
	
echo "Analysis completed." >> $LOG;
