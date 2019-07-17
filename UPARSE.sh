#!/bin/bash

set -o errexit

echo -e "Uparse pipeline started at $(date|awk '{print $4}')\nUsearch version 8.0.1623" >> $LOG;
		
		
#Lets call the function to check the dimension of multiplexed/multiplexed_linearized.fasta files and store it in a variable
check_dimension

# Set a variable to give a short name for the USEARCH binary
UPARSE=$(which usearch8.0.1623_i86linux32);
	
mkdir -p Multiplexed_files/Uparse_pipeline
# If the dimension is in megabyte or is smaller then 3 gigabyte we can use usearch for dereplication
if [[ -n "$megabyte" ]] || [[ -n "$gigabyte" ]] && [[ "$gigabyte" -le "3" ]]; then
		echo "-derep_fulllength $(date|awk '{print $4}')" >> $LOG;
		$UPARSE -derep_fulllength Multiplexed_files/multiplexed_linearized.fasta -fastaout Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated.fasta -sizeout -threads $NSLOTS;
		SING=$(grep singletons *.e$JOB_ID);
		echo $SING| awk '{print "singletons discarded "$9" ("$7" reads)"}'|sed 's/(//;s/)//' >> $LOG ;
else
		# I will use a piece of Vsearch just for demultiplexing.
		#grep -v "^>" Multiplexed_files/multiplexed_linearized.fasta | grep -v [^ACGTacgt] | sort -d | uniq -c | while read abundance sequence ; do hash=$(printf "${sequence}" | sha1sum); hash=${hash:0:40};printf ">%s;size=%d;\n%s\n" "${hash}" "${abundance}" "${sequence}"; done >> Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated.fasta
		echo "Uparse -derep_fulllength cannot handle files bigger than 3GB. I have to use a piece of Vsearch just for this task."
		echo "--derep_fulllength $(date|awk '{print $4}')" >> $LOG;
		VSEARCH=$(which vsearch113);
		$VSEARCH --derep_fulllength Multiplexed_files/multiplexed_linearized.fasta --output Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated.fasta --sizeout;				
fi ;
		
# Filter singletons
$UPARSE -sortbysize Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated.fasta -minsize 2 -fastaout Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2.fasta;

# cluster OTUs and de novo chimera removal
echo "cluster_otus (de novo chimera removal) $(date|awk '{print $4}')" >> $LOG;
$UPARSE -cluster_otus Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2.fasta -otus Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset.fasta;
CHIM_DE=$(awk '/chimeras/' *.e$JOB_ID);
echo "{$CHIM_DE}"|head -1|awk 'END {print "De novo chimera removal found : "$NF" Chimeras ("$(NF-2)" OTUs). "$(NF-4)" OTUs left."}'| sed 's/(//1;s/)//1'|sed 's/\r//g' >> $LOG;
		
# Reference chimera removal
echo "-uchime_ref (reference chimera removal) $(date|awk '{print $4}')" >> $LOG;
$UPARSE -uchime_ref Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset.fasta -db $CHIM -strand plus -nonchimeras Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras.fasta -threads $NSLOTS;
#CHIM_REF=$(awk '/chimeras/' *.e$JOB_ID);
#echo "Ref. based chimera check found :" >> $LOG;
CHIM_REF=$(grep -E chimeras\|found *.e$JOB_ID);
#I am not sure this will work any time
echo "{$CHIM_REF}"|grep -v +|head -2|tail -1|awk 'END {print "Reference based chimera removal found : "$NF" Chimeras ("$(NF-3)" OTUs)"}'| sed 's/(//1;s/)//1'|sed 's/\r//g' >> $LOG;
		
# Label OTUs using UPARSE python script
#python $( which fasta_number.py) Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras.fasta OTU_> Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta;

python $DIR/relabel_fasta/relabel_fasta.py Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras.fasta OTU_> Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta;

# Map reads (including singletons) back to OTUs directly when possible. If the file is too big the script will use an alternative strategy from Dr Ijaz tutorial
mkdir -p Uparse_OTU_tables;		
if [[ -n "$megabyte" ]] || [[ -n "$gigabyte" ]] && [[ "$gigabyte" -le "3" ]]; then
		echo "-usearch_global $(date|awk '{print $4}')" >> $LOG;
		$UPARSE -usearch_global Multiplexed_files/multiplexed_linearized.fasta -db Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -strand plus -id 0.97 -uc Uparse_OTU_tables/otu.map.uc -threads $NSLOTS;

else
		#The alternative strategy consist of splitting the file in 100 parts, perform the analysis on every part and merge the results
		mkdir -p Uparse_OTU_tables/split_files;
		mkdir -p Uparse_OTU_tables/uc_files;
		perl $(which fasta-splitter.pl) -n-parts-total 100 --out-dir Uparse_OTU_tables/split_files/ Multiplexed_files/multiplexed_linearized.fasta;
set +e		
for i in $(ls Uparse_OTU_tables/split_files/*.fasta 2> /dev/null); do 
set -e
				$UPARSE -usearch_global $i -db Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -strand plus -id 0.97 -uc Uparse_OTU_tables/uc_files/$(basename ${i}).map.uc 2> /dev/null; 
set +e
		done;
set -e
		cat Uparse_OTU_tables/uc_files/* > Uparse_OTU_tables/otu.map.uc ;
fi;
			
# Make OTU table. This is not the original uc2otutab.py. It has been modified to get a different output. 
python $DIR/uc2otutab/uc2otutab.py Uparse_OTU_tables/otu.map.uc > Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.txt;

# Convert to biom
echo "Create otu table $(date|awk '{print $4}')" >> $LOG;
if [ -f Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom ]; then
		rm Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom;
fi;
biom convert --table-type="OTU table" -i Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.txt -o Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom --to-json;

# Assign taxonomy
#if [[ -z $SILVA ]]; then
#	echo "parallel_assign_taxonomy_rdp.py $(date|awk '{print $4}')" >> $LOG;
#	parallel_assign_taxonomy_rdp.py --rdp_max_memory 12000 -t $TAX -r $REF -i Multiplexed_files/Vsearch_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Vsearch_OTU_tables/assigned_taxonomy -O $NSLOTS;
#else		
#	echo "parallel_assign_taxonomy_uclust.py $(date|awk '{print $4}')" >> $LOG;
#	parallel_assign_taxonomy_uclust.py -t $TAX -r $REF -i Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Uparse_OTU_tables/assigned_taxonomy -O $NSLOTS;
#fi;
		
# Assign taxonomy
echo -e "assign_taxonomy_rdp.py $(date|awk '{print $4}')\nRDP_Classifier version 2.2" >> $LOG;
assign_taxonomy.py -m rdp --rdp_max_memory $MEMORY -t $TAX -r $REF -i Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Uparse_OTU_tables/assigned_taxonomy;
		
#add taxonomy to BIOM table
if [ -f Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom ]; then
		rm Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom;
fi;
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp Uparse_OTU_tables/assigned_taxonomy/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_tax_assignments.txt -i Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTU_table.biom -o Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom ;

#Filter out low abundance OTUs
filter_otus_from_otu_table.py -i Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom -o Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table_low_filtered.biom --min_count_fraction 0.00005 ;
		
#Align sequences using muscle. If you choose to use greengenes the variable $CORE is the same of $ALIGNED
echo "Create tree $(date|awk '{print $4}')" >> $LOG;
align_seqs.py -i Multiplexed_files/Uparse_pipeline/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs.fasta -o Uparse_OTU_tables/alignment/otus_aligned -t $CORE ;
		
# Filter the alignment for gaps
filter_alignment.py -i Uparse_OTU_tables/alignment/otus_aligned/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_aligned.fasta -o Uparse_OTU_tables/alignment ;
	
#Make tree. Required for the next step
make_phylogeny.py -i Uparse_OTU_tables/alignment/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_OTUs_aligned_pfiltered.fasta -o Uparse_OTU_tables/otus.tre ;
		
#Rename the otu table
mv Uparse_OTU_tables/multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table_low_filtered.biom Uparse_OTU_tables/otu_table.biom ;
		
#OTU table and tree UPARSE. The script will use these parameters during the last step
export BIOM="Uparse_OTU_tables/otu_table.biom" ;
export TREE="Uparse_OTU_tables/otus.tre" ;
		
echo "UPARSE pipeline finished at $(date|awk '{print $4}')" >> $LOG;
		
		
#Lets create the directory for the outputs in the third step
mkdir -p RESULTS/ ;
				
if [[ -z $SILVA ]]; then			
		mkdir -p RESULTS/Uparse_gg ;
		export RESULTS_PATH=RESULTS/Uparse_gg ;			
else			
		mkdir -p RESULTS/Uparse_silva ;
		export RESULTS_PATH=RESULTS/Uparse_silva ;			
fi;
