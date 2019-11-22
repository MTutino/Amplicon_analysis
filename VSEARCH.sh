#!/bin/bash

set -o errexit

echo -e "Vsearch pipeline started at $(date|awk '{print $4}')\nVsearch version 2.10.4" >> $LOG;

#**************JUST A TEST*********
#CORE="$DIR/Silva/Silva119_release_aligned_rep_files/90_16S_only/Silva_119_rep_set90_aligned_16S_only_gfiltered/Silva_119_rep_set90_aligned_16S_only_pfiltered.fasta"
#**********************************
		
# Set a variable to give a short name for the VSEARCH binary
#VSEARCH=$(which vsearch113);
VSEARCH=$(which vsearch);
which vsearch >> $LOG
	
mkdir -p Multiplexed_files/Vsearch_pipeline ;
# Dereplication
echo "--derep_fulllength (dereplicate and discard singletons) $(date|awk '{print $4}')" >> $LOG;
$VSEARCH --derep_fulllength Multiplexed_files/multiplexed_linearized.fasta --output Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2.fasta --sizeout --minuniquesize 2;
SINGLETONS=$(grep "uniques written" *.e$JOB_ID);
echo $SINGLETONS|awk '{print "Singletons discarded "$7" ("$4" reads)"}'|sed 's/(//;s/)//' >> $LOG;
				
# De novo chimera removal
echo "--uchime_denovo $(date|awk '{print $4}')" >> $LOG;
$VSEARCH --uchime_denovo Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2.fasta --nonchimeras Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_DNchimerafiltered.fasta;

# cluster OTUs
echo "--cluster_fast $(date|awk '{print $4}')" >> $LOG;
$VSEARCH --cluster_size Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_DNchimerafiltered.fasta --centroids Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_DNchimerafiltered.fasta --id 0.97 --threads $NSLOTS;
				
# Reference chimera removal
echo -e "--uchime_ref $(date|awk '{print $4}')\nReference: RDPClassifier_16S_trainsetNo14" >> $LOG;
$VSEARCH --uchime_ref Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_DNchimerafiltered.fasta --db $CHIM  --nonchimeras Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras.fasta;				
CHIMERAS=$(grep -E Found\|suspicious *.e$JOB_ID);
echo $CHIMERAS|awk '{print "De novo chimera removal found "$3" chimeras ["$2" seqs], "$6" non chimeras ["$5" seqs], suspicious candidates "$10" ["$9" seqs]\nClosed-reference chimera removal found "$18" chimeras ["$17" OTUs], "$21" non chimeras ["$20" OTUs], suspicious candidates "$25" ["$24" OTUs]"}'|sed 's/(//g;s/)//g' >> $LOG;
#echo $CHIMERAS|awk '{print "De novo chimera removal found "$9" chimeras ["$8" seqs], "$12" non chimeras ["$11" seqs], suspicious candidates "$16" ["$15" seqs]\nClosed-reference chimera removal found "$24" chimeras ["$23" OTUs], "$27" non chimeras ["$26" OTUs], suspicious candidates "$31" ["$30" OTUs]"}'|sed 's/(//g;s/)//g' >> $LOG;	
#Before adding the metatable check
#CHIMERAS=$(grep -E Found\|suspicious *.e$JOB_ID);
#echo $CHIMERAS|awk '{print "De novo chimera removal found "$3" chimeras, "$6" non chimeras, suspicious candidates "$10"\nClosed-reference chimera removal found "$18" chimeras, "$21" non chimeras, suspicious candidates "$25}'|sed 's/(//g;s/)//g' >> $LOG;
				
# Label OTUs 
$VSEARCH --sortbysize Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras.fasta --relabel OTU_ --output Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta;

# Map reads (including singletons) back to OTUs 
# vsearch version > 2.2 has option --biomout to output directly a biom table
mkdir -p Vsearch_OTU_tables
$VSEARCH --usearch_global Multiplexed_files/multiplexed_linearized.fasta --db Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta --strand plus --id 0.97 --uc Vsearch_OTU_tables/otu.map.uc;

# Make OTU table
echo "Create otu table $(date|awk '{print $4}')" >> $LOG;
python $DIR/uc2otutab/uc2otutab.py Vsearch_OTU_tables/otu.map.uc > Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.txt;

# Convert to biom
if [ -f Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom ]; then
		rm -f Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom;
fi;
biom convert --table-type="OTU table" -i Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.txt -o Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom --to-json;


#This is a try with a reference trimmed to the V4 region only.
#REF="$DIR/Silva/Silva119_release/rep_set/97/97_trimmed_2.fasta"



# Assign taxonomy
echo -e "assign_taxonomy_rdp.py $(date|awk '{print $4}')\nRDP_Classifier version 2.2" >> $LOG;
assign_taxonomy.py -m rdp --rdp_max_memory $MEMORY -t $TAX -r $REF -i Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Vsearch_OTU_tables/assigned_taxonomy ;
				
#add taxonomy to BIOM table. Remove a previous biom table first because "biom" do not overwrite the file and it would return an error.
if [ -f Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom ]; then
		rm -f Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom;
fi;
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp Vsearch_OTU_tables/assigned_taxonomy/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_tax_assignments.txt -i Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom -o Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom ;

#Filter out low abundance OTUs
filter_otus_from_otu_table.py -i Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom -o Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table_low_filtered.biom --min_count_fraction 0.00005 ;
				
echo "Create tree $(date|awk '{print $4}')" >> $LOG;
#Align sequences using muscle. If you chose to use greengenes the variable $CORE is the same of $ALIGNED
align_seqs.py -i Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Vsearch_OTU_tables/alignment/otus_aligned -t $CORE ;
		
# Filter the alignment for gaps
filter_alignment.py -i Vsearch_OTU_tables/alignment/otus_aligned/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_aligned.fasta -o Vsearch_OTU_tables/alignment ;
	
#Make tree. Required for the next step
make_phylogeny.py -i Vsearch_OTU_tables/alignment/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_aligned_pfiltered.fasta -o Vsearch_OTU_tables/otus.tre ;
		
#Rename the otu table
mv Vsearch_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table_low_filtered.biom Vsearch_OTU_tables/otu_table.biom ;
		
#OTU table and Tree VSEARCH. The script will use these parameters during the last step.
export BIOM="Vsearch_OTU_tables/otu_table.biom" ;
export TREE="Vsearch_OTU_tables/otus.tre" ;
		
		
echo "Vsearch pipeline finished at $(date|awk '{print $4}')" >> $LOG;
		
		
#Lets create the directory for the outputs in the third step
mkdir -p RESULTS/ ;
				
if [[ -z $SILVA ]]; then
		mkdir -p RESULTS/Vsearch_gg ;
		export RESULTS_PATH="RESULTS/Vsearch_gg" ;			
else			
		mkdir -p RESULTS/Vsearch_silva ;
		export RESULTS_PATH="RESULTS/Vsearch_silva" ;		
fi;
