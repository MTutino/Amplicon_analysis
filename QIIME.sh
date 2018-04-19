#!/bin/bash

set -o errexit

echo -e "QIIME pipeline started at $(date|awk '{print $4}') \nQIIME version 1.8.0" >> $LOG;
		
# QIIME pipeline with Usearch61

#Lets define the configuration of QIIME		
if [[ -z $SILVA ]]; then
		awk -v k=$REF_DATA_PATH '{if (NR==12) $0="pynast_template_alignment_fp "k"/gg_13_8_otus/rep_set_aligned/97_otus.fasta"} {if (NR==23) $0="assign_taxonomy_reference_seqs_fp "k"/gg_13_8_otus/rep_set/97_otus.fasta"} {if (NR==24) $0="assign_taxonomy_id_to_taxonomy_fp "k"/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt"} {print}' $DIR/Config_QIIME/.qiime_config > ~/.qiime_config ;
else	
		awk -v k=$REF_DATA_PATH '{if (NR==12) $0="pynast_template_alignment_fp "k"/Silva/Silva119_release/core_alignment/core_Silva119_alignment.fna"} {if (NR==23) $0="assign_taxonomy_reference_seqs_fp "k"/Silva/Silva119_release/rep_set/97/Silva_119_rep_set97.fna"} {if (NR==24) $0="assign_taxonomy_id_to_taxonomy_fp "k"/Silva/Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt"} {print}' $DIR/Config_QIIME/.qiime_config > ~/.qiime_config;
fi;
				
# Fix the path to the qiime_scripts_dir in the QIIME config file
sed -i 's,\(qiime_scripts_dir\t\).*,\1'$(dirname $(which print_qiime_config.py))',g' ~/.qiime_config

#Lets call the function to check the dimension of multiplexed/multiplexed_linearized.fasta files and store it in a variable
check_dimension


# If the dimension is in megabyte or is smaller then 3 gigabyte we can use usearch61...
if [[ -n "$megabyte" ]] || [[ -n "$gigabyte" ]] && [[ "$gigabyte" -le "3" ]]; then
		echo "I can use Usearch version 6.1.544" >> $LOG;
		
		# QIIME pipeline open-reference
		
		# Identify chimeric seqs with Usearch61 (open reference)
		echo "identify_chimeric_seqs.py $(date|awk '{print $4}')" >> $LOG;
		identify_chimeric_seqs.py -m usearch61 -i Multiplexed_files/multiplexed_linearized.fasta -r $REF -o  usearch61_chimera_checking/ ;
		DENOVO_CHIM=$(grep denovo_chimeras usearch61_chimera_checking/identify_chimeric_seqs.log|awk '{print $2}');
		DENOVO_NON_CHIM=$(grep denovo_non_chimeras usearch61_chimera_checking/identify_chimeric_seqs.log|awk '{print $2}');
		REF_CHIM=$(grep ref_chimeras usearch61_chimera_checking/identify_chimeric_seqs.log|awk '{print $2}');
		REF_NON_CHIM=$(grep ref_non_chimeras usearch61_chimera_checking/identify_chimeric_seqs.log|awk '{print $2}'); 
		grep -v uchime_denovo  usearch61_chimera_checking/multiplexed_linearized.fasta_chimeras_denovo.log|awk '/chimeras/ {print "Uchime de-novo found "$(NF - 2)" chimeric OTUs "$NF"..."}' >> $LOG;
		echo  "..."$DENOVO_CHIM" sequences on "$DENOVO_NON_CHIM >> $LOG;
		grep -v uchime_ref  usearch61_chimera_checking/multiplexed_linearized.fasta_chimeras_ref.log|awk '/chimeras/ {print "Uchime ref. found "$(NF - 2)" chimeric OTUs "$NF}' >> $LOG;
		echo "..."$REF_CHIM" sequences on "$REF_NON_CHIM >> $LOG ;
		# Filter out chimeric seqs from Fasta file
		filter_fasta.py -f Multiplexed_files/multiplexed_linearized.fasta -o usearch61_chimera_checking/seqs_chimeras_filtered.fasta -s usearch61_chimera_checking/chimeras.txt -n ;
		
		# Create a parameter file for "pick_otus.py"	
		rm -f QIIME_OTU_tables/open_ref/ucr_open_ref_params.txt;
		rm -f ucr_open_ref_params.txt;
		echo -e pick_otus:otu_picking_method usearch61"\n"pick_otus:enable_rev_strand_match True"\n"pick_otus:max_rejects 64"\n"assign_taxonomy:assignment_method rdp"\n"assign_taxonomy:rdp_max_memory $MEMORY >> ucr_open_ref_params.txt;
		
		# Pick otus using an open reference method (Usearch61)
		echo "pick_open_reference_otus.py $(date|awk '{print $4}')" >> $LOG;
		rm -f QIIME_OTU_tables/open_ref/otu_table_mc2_w_tax.biom ;
		rm -f QIIME_OTU_tables/open_ref/otu_table_mc2_w_tax_no_pynast_failures.biom ;
		pick_open_reference_otus.py -i usearch61_chimera_checking/seqs_chimeras_filtered.fasta -o open_ref/ -r $REF -p ucr_open_ref_params.txt -a -O $NSLOTS -f ;

		#Filter out low abundance OTUs
		filter_otus_from_otu_table.py -i open_ref/otu_table_mc2_w_tax_no_pynast_failures.biom -o open_ref/otu_table_no_chimeras_no_singletons_no_low_abund_OTUs.biom --min_count_fraction 0.00005 ;
				
		mkdir -p QIIME_OTU_tables ;
		mkdir -p QIIME_OTU_tables/open_ref ;
		rsync -a open_ref/ QIIME_OTU_tables/open_ref/ ;
		rm -r open_ref/ ;
		mv ucr_open_ref_params.txt QIIME_OTU_tables/open_ref ;
				
		#Rename the otu table
		mv QIIME_OTU_tables/open_ref/otu_table_no_chimeras_no_singletons_no_low_abund_OTUs.biom QIIME_OTU_tables/open_ref/otu_table.biom ;
				
		#OTU table and tree QIIME open reference. The script will use these parameters during the last step
		export BIOM="QIIME_OTU_tables/open_ref/otu_table.biom" ;
		export TREE="QIIME_OTU_tables/open_ref/rep_set.tre" ;
								
else
		# ... otherwise we have to use an alternative method
		# QIIME pipeline closed-reference
		echo "I am sorry, I cannot use Usearch61, multiplexed_linearized.fasta is too big. I will use a closed-reference method." >> $LOG;

		if [[ ! -z $SILVA ]]; then
				echo "ChimeraSlayer does not support the latest version of Silva db yet. For chimera detection Silva version 108 will be used." >> $LOG;
				ALIGNED="$REF_DATA_PATH/Silva/Silva_108_core_alignment/Silva_108_core_aligned_seqs.fasta"
		fi;

		# Create a parameter file for "pick_otus.py".rm if it already exists
		rm -f QIIME_OTU_tables/closed_ref/ucr_closed_ref_params.txt;
		echo -e "\n"pick_otus:enable_rev_strand_match True"\n"pick_otus:max_accepts 1"\n"pick_otus:max_rejects 8"\n"pick_otus:stepwords 8"\n"pick_otus:word_length 8 >> ucr_closed_ref_params.txt;
			
		# Pick otus using a closed-reference method (Uclust) 
		echo "pick_closed_reference_otus.py $(date|awk '{print $4}')" >> $LOG;
		pick_closed_reference_otus.py -i Multiplexed_files/multiplexed_linearized.fasta -r $REF -o closed_ref/ -t $TAX -p ucr_closed_ref_params.txt -a -O $NSLOTS -f;
		mv closed_ref/otu_table.biom closed_ref/otu_table_with_chimeras.biom
			
		#Pick rep set
		echo "pick_rep_set.py with Uclust $(date|awk '{print $4}')" >> $LOG;
		pick_rep_set.py -i closed_ref/uclust_ref_picked_otus/multiplexed_linearized_otus.txt -r $REF -f Multiplexed_files/multiplexed_linearized.fasta -o Multiplexed_files/multiplexed_linearized.fasta_rep_set.fasta;
			
		#Align sequences with Pynast
		echo "align_seqs.py $(date|awk '{print $4}')" >> $LOG;
		align_seqs.py -i Multiplexed_files/multiplexed_linearized.fasta_rep_set.fasta -t $ALIGNED ;
			
		#Identify chimeras with Chimeraslayer
		echo "identify_chimeric_seqs.py  with ChimeraSlayer (r20110519) $(date|awk '{print $4}')" >> $LOG;
		mkdir -p Chimeraslayer_checking/ ;
		identify_chimeric_seqs.py -m ChimeraSlayer -i pynast_aligned/multiplexed_linearized.fasta_rep_set_aligned.fasta -a $ALIGNED -o Chimeraslayer_checking/multiplexed_linearized.fasta_rep_set_aligned_chimeric.txt ;
				
		#Singleton filtering
		filter_otus_from_otu_table.py -i closed_ref/otu_table_with_chimeras.biom -o closed_ref/otu_table_with_chimeras_no_singletons.biom -n 2 ;
			
		#Filter out chimeras from the otu table
		filter_otus_from_otu_table.py -i closed_ref/otu_table_with_chimeras_no_singletons.biom -o closed_ref/otu_table_no_chimeras_no_singletons.biom -e Chimeraslayer_checking/multiplexed_linearized.fasta_rep_set_aligned_chimeric.txt;
		
		#Filter out low abundance OTUs
		filter_otus_from_otu_table.py -i closed_ref/otu_table_no_chimeras_no_singletons.biom -o closed_ref/otu_table_no_chimeras_no_singletons_no_low_abund_OTUs.biom --min_count_fraction 0.00005 ;
				
		mkdir -p QIIME_OTU_tables ;
		mkdir -p QIIME_OTU_tables/closed_ref ;
		rsync -a closed_ref/ QIIME_OTU_tables/closed_ref ;
		rm -r closed_ref/ ;
		mv ucr_closed_ref_params.txt QIIME_OTU_tables/closed_ref ; 
				
		#Rename the otu table
		mv QIIME_OTU_tables/closed_ref/otu_table_no_chimeras_no_singletons_no_low_abund_OTUs.biom QIIME_OTU_tables/closed_ref/otu_table.biom ;
				
		#OTU table end tree QIIME closed reference. The script will use these parameters during the last step
		export BIOM="QIIME_OTU_tables/closed_ref/otu_table.biom" ;
	
fi;
		
echo "QIIME pipeline finished at $(date|awk '{print $4}')" >> $LOG;
		
		
#Lets create the directory for the outputs in the third step
		
mkdir -p RESULTS/ ;
				
if [[ -n $SILVA ]]; then			
	mkdir -p RESULTS/QIIME_silva ;
	export RESULTS_PATH=RESULTS/QIIME_silva ;
	
elif [[ -n $HOMD ]];then
	mkdir -p RESULTS/QIIME_homd ;
	export RESULTS_PATH=RESULTS/QIIME_homd ;

else			
	mkdir -p RESULTS/QIIME_gg ;
	export RESULTS_PATH=RESULTS/QIIME_gg ;				
fi;
