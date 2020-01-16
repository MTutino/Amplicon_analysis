**This script has been designed for the analysis of 16S rRNA data from Illumina Miseq (Casava >= 1.8) paired-end reads** 

It is divided in 3 steps:

	1) Quality control with Cutadapt to check for and trim PCR primers, Sickle for quality trimming, SPADes for illumina Miseq errors correction and Pandaseq to merge paired-end reads
	    (Based on the Dr.UZ Ijaz’s command line tutorial “Illumina Amplicons OTU Construction with Noise Removal” found at www.tinyurl.com/JCBioinformatics)
	    (Schirmer and Ijaz et al.,"Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform." ,Nucleic Acids Res., 2015)
	    
	2) Remove singletons and chimeras, and build an OTU table
		REMOVED: QIIME pipeline: to remove singletons and chimeras and build otu table and phylogenetic tree (open-reference method with usearch61 if possible (it depends on the files' size), otherwise closed-reference method with #		uclust)
	
		UPARSE pipeline: free version limited to 4Gb RAM
	
		Vsearch pipeline: A freely available programme almost identical to Uparse but not limited to 4Gb RAM
		
		DADA2: An R package to infer exact amplicon sequence variants (ASVs) 
	
	3)  QIIME to perform beta and alpha diversity analysis




The script requires the files "Final_name.txt" and "Metatable.txt". "Categories.txt" is optional.
Pay attention to the name of these files, they have to match!

You can find an example of these files below.

Create a folder (the name is not important) with raw files to analise inside (fastq or fastq.gz files), Final_name.txt, Metatable.txt and Categories.txt. 
You will have to run the programme (or the jobscript if you are using it in CSF) inside that folder.


#REQUIRED FILES

**Final_name.txt:**

Final_name.txt is a text tab-delimited file of two columns without any header. 
The first column is a list of your files' name, you can get the list this way

	for i in $(ls *fastq*);do echo $(basename ${i});done

The second column is a list of names decided by the user (final name). The final name will be the one you will see in the final plots, decide a meaningful name. Every final name has to be unique and specific for 
a sample. Forward (R1) and reverse (R2) reads from the same sample have to be linked to the same final name.  
For the final name the only accepted special characters are underscore "_", dot "." and dash "-".

**_You can create this file in Excel but be careful to save it as Text(Tab delimited)_**
**_After importing the file from Windows run the command "dos2unix Final_name.txt" to convert it in unix format_**

**EXAMPLE:**

	Mock_S4_L001_R1_001.fastq.gz	Mock_RUN1
	Mock_S4_L001_R2_001.fastq.gz	Mock_RUN1
	Mock_S5_L001_R1_001.fastq.gz 	Mock_RUN2
	Mock_S5_L001_R2_001.fastq.gz	Mock_RUN2
	Mock_S6_L001_R1_001.fastq.gz  	Mock_RUN3
	Mock_S6_L001_R2_001.fastq.gz	Mock_RUN3


**Metatable.txt:**

This file is important for the third step. You have to use a specific format. 
You can find a description of it on QIIME website for more information (http://qiime.org/documentation/file_formats.html).
The column Run, defining the sequencing run, must be included when using DADA2. If the samples were sequenced in the same run, simply use the same run identifier.

**EXAMPLE:**

	#SampleID	BarcodeSequence	LinkerPrimerSequence	Run	Description
	Mock-RUN1	TAAGGCGAGCGTAAGA		1	Control
	Mock-RUN2	CGTACTAGGCGTAAGA		2	Control
	Mock-RUN3	AGGCAGAAGCGTAAGA		3	Control

The column "LinkerPrimerSequence" is empty but it cannot be deleted.
The header is very important. "#SampleID", "Barcode", "LinkerPrimerSequence", "Run" and "Description" are mandatory. 
Between "LinkerPrimerSequence" and "Description" you can add as many columns as you want.
For every column a PCoA plot will be created during the third step.

**_You can create this file in Excel but be careful to save it as Text(Tab delimited)_**
**_After importing the file from Windows run the command "dos2unix Metatable.txt" to convert it in unix format_**

**_During the analysis the Metatable.txt will be checked to be sure that the metatable has the correct format_**
**_If necessary this will be modified. You can find the new corrected metatable file in the folder "Metatable_log/Metatable_mod.txt"_** 
**_If you are going to use the metatable file for other statistical analyses remind to use the modified one, otherwise the samples name will not match!_**

**OPTIONAL FILE:**

**Categories.txt:**

This file is required if you want to get box plots for the alpha diversity comparison. You will get a box plot for every category you write in the file with different metrics.
It is just a list (without header and IN ONE COLUMN) of categories present in the Metatable.txt file.
THE NAMES YOU ARE USING HAVE TO BE THE SAME YOU USED IN THE METATABLE.TXT

**_You can create this file in Excel but be careful to save it as Text(Tab delimited)_**
**_After importing the file from Windows run the command "dos2unix Categories.txt" to convert it in unix format_**


**EXAMPLE** (these categories are useless):

	BarcodeSequence
	Run
	LinkerPrimerSequence
	Description


Now that you created all the required files, and copied the raw files, let's see how to use the programme

#MODULES
If you are using the programme on CSF load the following modules (you can copy and paste):

**Step1:**

	#(Cutadapt is included in anaconda)
	module load apps/anaconda/2.5.0/bin
	module load apps/gcc/sickle/1.33
	module load apps/gcc/bioawk/27-08-2013
	module load apps/pandaseq/2.8.1/gcc-4.8.5
	module load apps/binapps/spades/3.10.1
	module load apps/fastqc/0.11.3/noarch

**Step2/3:**

	module load apps/bioinf
	module load apps/qiime/1.9.1/python-2.7.8+numpy-1.9.2+scipy-0.17.0+pycogent-1.5.3+pandas-0.17.0+biomformat-2.1.4+matplotlib-1.4.3+ghc-7.8.2+pynast-1.2.2
	module load apps/binapps/vsearch/2.10.4
	module load apps/binapps/fasta-splitter/0.2.6
	module load apps/rdpclassifier/2.2/noarch
	module load compilers/gcc/8.2.0
	module load apps/binapps/blast/legacy/2.2.26
	module load apps/binapps/usearch/8.0.1623
	module load apps/R/3.5.2/gcc-4.8.5+lapack-3.5.0+blas-3.6.0


On the other hand, if you are using the programme in interactive mode be sure you have installed all the required programmes and that they are installed in your bin folder. 


#OPTIONS
Best options for the analysis of V3-V4 hypervariable regions [Default options]: -q 20 -l 10 -o 10 -L 380

	-h      Show this message
   	-g      Forward PCR primer, without any barcode/adapter (Cutadapt) ***REQUIRED IF USING CUTADAPT. IF IT IS NOT PASSED, CUTADAPT WILL BE DISABLED***
   	-G      Reverse PCR primer, without any barcode/adapter (Cutadapt) ***REQUIRED IF USING CUTADAPT.IF IT IS NOT PASSED, CUTADAPT WILL BE DISABLED***
   	-q      Phred score threshold below which the read will be trimmed [default 20] (Sickle) ***OPTIONAL*** 
   	-l (Lowercase "L")	   Length of the sliding Window in bp [default 10] (Sickle) ***OPTIONAL***
   	-O (Uppercase "O")     Minimum overlap in bp between forward and reverse reads [default 10] (Pandaseq) ***OPTIONAL***
   	-L      Minimum length in bp for a sequence to be kept after overlapping [default 380] (Pandaseq) ***OPTIONAL***
   	-1 (One)     Use this option "-1 suppress" to skip the Quality Control step
   	-P	   Use this option to decide which pipeline you want to use, UPARSE, Vsearch or DADA2. UPARSE="-P uparse". Vsearch="-P vsearch". DADA2="-P DADA2"  ***REQUIRED***
   	-S	   The default reference database is GreenGenes. Use this option without any argument if you want to use Silva. To use Silva you need at least 22 Gb of RAM.

 ***_To run only the third step_***
 
   	-3	   Pass this flag without any argument to run only the third step
   	-i	   BIOM file
   	-o (Lowercase "O")   Path to the directory where you want to store the result
   	-m	   Metatable file. The names in the metatable has to match with those used to create the BIOM file
   	-t	   Tree file. It has to be same as used in the pipeline (cannot use SILVA to analyse the data and GreenGenes for the plots)
   
**_IF YOU WANT TO USE SILVA DO NOT RUN THE SCRIPT WITH LESS THAN 5 CORES (at least 22 gb of RAM)_**
**_QIIME with a large dataset can be very slow. It can take few days to finish the analysis. I advice you to use either Uparse or Vsearch_**

 ***_Advanced options_***

Only use these if you know what you're doing and exercise caution!

   	-r	   Path to the directory with the reference databases, if not same as the script directory. If used then the Greengenes and Silva database directories must be subdirectories of this path. ***OPTIONAL***

#EXAMPLE USAGE IN THE INTERACTIVE MODE:

1)The following command will use default options for QC, the Vsearch pipeline and greengenes database

	Full_Path_to_the_script/amplicon_analysis_pipeline.sh -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -P vsearch 

2)The following command will NOT use Cutadapt (e.g. to analyse samples amplified with different primers together. Not raccomanded), will use the uparse pipeline and silva database

	Full_Path_to_the_script/Amplicon_analysis_pipeline.sh -P uparse -S

3)The following command will suppress the first step and the analysis will start directly from the second step.
Useful if you have already analysed these data, for example with Vsearch, and you want to try a different pipeline or if you used GreenGenes and you want to try Silva

	Path_to_the_script/Amplicon_analysis_pipeline.sh -1 suppress -P uparse

4)The following command will run only the third step
Useful if you want to produce new plots after filtering out some samples from the OTU table.

	Path_to_the_script/Amplicon_analysis_pipeline.sh -3 -i Path/Vsearch_OTU_tables/otu_table_without_mocks.biom -o Path/RESULTS_2/ -m Path/Metatable_log/Metatable_mod.txt -t Path/Vsearch_OTU_tables/otus.tre

#EXAMPLE USAGE IN CSF:

**To run this programme in CSF you have to write a jobscript. To do that you can use gedit.
Type gedit on the command line, this will open the text editor, then copy and paste the following lines
Change the path and the flags as explained above**

	#!/bin/bash -x
	#$ -S /bin/bash   # Inform SGE we are using the bash shell
	#$ -cwd           # Job will run in the current directory (where you ran qsub)
	#$ -V             # Inherit current environment (e.g., any loaded modulefiles)
                  # ... important so that commands can be found when jobs run.
	#$ -pe smp.pe 5   # Define how many cores you want to use. REMEMBER NOT TO USE LESS THEN FIVE CORES FOR SILVA 

	** #The following command will use the vsearch pipeline and greengenes as database**
	Path_to_the_script/mplicon_analysis_pipeline.sh -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -P vsearch


**IF YOU WANT TO REPEAT THE ANALYSIS FILTERING OUT SOME SAMPLES YOU WILL HAVE TO MOVE,RENAME OR DELETE THE "MULTIPLEXED_FILES" FOLDER, THOSE FILES ARE NEVER OVERWRITTEN. 
ON THE OTHER SIDE, ALL THE OTHER FILES WILL BE OVERWRITTEN SO, IF YOU RUN THE PIPELINE WITH THE SAME PARAMETERS, MOVE OR RENAME THE FILES YOU WOULD LIKE TO KEEP BEFORE RUNNING THE PROGRAMME.**

If the programme ran successfully you will get the following output


#OUTPUT:

**"QUALITY_CONTROL"**: It contains information about length distribution for raw and trimmed data plus the number of reads for every piece of the quality control step.


**"DADA2_OTU_tables"** : 
		1) Error_rate_plots # Folder containing the plots of learned error rates
		2) DADA2_tax_OTU_table.biom #Otu table not filtered for low abundance OTUs
		3) DADA2_tax_OTU_table_low_filtered.biom #Otu table filtered for low abundance OTUs. filter_otus_from_otu_table.py -i INPUT -o OUTPUT --min_count_fraction 0.00005
		4) otus.tre #The tree file


**"Vsearch_OTU_tables"**, **"Uparse_OTU_tables"**: These folders have the same files inside.

		1) multiplexed_linearized_dereplicated_mc2_repset_nonchimeras_tax_OTU_table.biom. #Otu table not filtered for low abundance OTUs
		2) otu_table.biom. Otu table filtered for low abundance OTUs. #filter_otus_from_otu_table.py -i INPUT -o OUTPUT --min_count_fraction 0.00005
		3) otus.tre. #The tree file





**"RESULTS"**: Inside this folder you will find a subfolder for every pipeline/reference DB you used (e.g. QIIME_silva/  Uparse_silva/  Vsearch_gg/  Vsearch_silva/) but the content will be the same.

	-table_summary.txt:
		A file with the sequences count per sample and the total number of OTUs	
	
	-OTUs_count.txt:
		A file with the number of unique observations per sample (rather than the total count of observations per sample)

	-taxa:
		An otu table for every level (Phylum, genus, species...)

	-Rarefied_otu_tables:
		Otu tables at different levels of rarefaction

	-phylum_genus_charts:
		bar_charts.html. An HTML file to to visualise the relative abundance of taxa in your samples. Only at a Phylum, Genus and Species level
		pie_charts.html. As the previous one.

	-Alpha_diversity:
		-collated_alpha: Each file represents one diversity metric, each row in a file represents one (rarefied) otu table and each column in a file represents one sample
		-Alpha_rarefactions: Alpha diversity on multiple OTU tables (Rarefied_otu_tables/) with different metrics (observed_species,chao1,shannon,simpson,fisher_alpha)
		-Alpha_diversity_boxplot: This folder is generated only if the file "Categories.txt" is present inside the working directory.
		 Inside this folder you will find a subfolder for every metric (observed_species,chao1,shannon,simpson,fisher_alpha) and inside those subfolders a box plot for 
		 every category you wrote in the file "Categories.txt".
		-rarefaction_curves: it containes an html file of rarefaction plots
		
	-beta_div_even:
		2-D PCoA plots from beta diversity analysis at the minimum depth with weighted and unweighted unifrac distance metrics 

	-Jackknifed_betadiversity:
		To directly measure the robustness of individual UPGMA clusters and clusters in PCoA plots with both unweighted_unifrac/ and weighted_unifrac/
		Inside the folders unweighted_unifrac/ and weighted_unifrac/ you will find a bootstrapped tree "jackknife_unweighted_unifrac.pdf"/"jackknife_weighted_unifrac.pdf"
		and a folder called emperor_pcoa_plots/ with 3-D PCOA plots. 

