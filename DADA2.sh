#!/bin/bash

set -o errexit

echo -e "DADA2 pipeline started at $(date|awk '{print $4}')\nDADA2 version 3.8" >> $LOG

mkdir -p DADA2_OTU_tables
mkdir -p DADA2_OTU_tables/Error_rate_plots


Rscript $DIR/DADA2.R

# Convert to biom
if [ -f DADA2_OTU_tables/DADA2_OTU_table.biom ]; then
		rm -f DADA2_OTU_tables/DADA2_OTU_table.biom;
fi;
biom convert --table-type="OTU table" -i DADA2_OTU_tables/seq_table.txt -o DADA2_OTU_tables/DADA2_OTU_table.biom --to-json

#add taxonomy to BIOM table. Remove a previous biom table first because "biom" do not overwrite the file and it would return an error.
if [ -f DADA2_OTU_tables/DADA2_tax_OTU_table.biom ]; then
		rm -f DADA2_OTU_tables/DADA2_tax_OTU_table.biom;
fi;
biom add-metadata --sc-separated taxonomy --observation-header OTUID,taxonomy --observation-metadata-fp DADA2_OTU_tables/seq_Taxonomy.txt -i DADA2_OTU_tables/DADA2_OTU_table.biom -o DADA2_OTU_tables/DADA2_tax_OTU_table.biom ;


#Filter out low abundance OTUs
filter_otus_from_otu_table.py -i DADA2_OTU_tables/DADA2_tax_OTU_table.biom -o DADA2_OTU_tables/DADA2_tax_OTU_table_low_filtered.biom --min_count_fraction 0.00005 ;
	
echo "Create tree $(date|awk '{print $4}')" >> $LOG;
#Align sequences using muscle. If you chose to use greengenes the variable $CORE is the same of $ALIGNED
align_seqs.py -i DADA2_OTU_tables/seqs.fa -o DADA2_OTU_tables/otus_aligned -t $CORE ;
# For testing 
#align_seqs.py -i DADA2_OTU_tables/seqs.fa -o DADA2_OTU_tables/otus_aligned -t /mnt/jw01-aruk-home01/projects/psa_microbiome/common_files/scripts/Amplicon_analysis_pipeline_Silva123/Silva/SILVA123_QIIME_release/core_alignment/core_alignment_SILVA123.fasta
	
# Filter the alignment for gaps
filter_alignment.py -i DADA2_OTU_tables/otus_aligned/seqs_aligned.fasta -o DADA2_OTU_tables/alignment ;

#Make tree. Required for the next step
make_phylogeny.py -i DADA2_OTU_tables/alignment/seqs_aligned_pfiltered.fasta -o DADA2_OTU_tables/otus.tre ;

#OTU table and Tree DADA2. The script will use these parameters during the last step.
export BIOM="DADA2_OTU_tables/DADA2_tax_OTU_table_low_filtered.biom" ;
export TREE="DADA2_OTU_tables/otus.tre" ;	
		
echo "DADA2 pipeline finished at $(date|awk '{print $4}')" >> $LOG;
			
#Lets create the directory for the outputs in the third step
mkdir -p RESULTS/ ;
mkdir -p RESULTS/DADA2_silva ;
export RESULTS_PATH="RESULTS/DADA2_silva" ;
