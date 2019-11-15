###############	LOAD PACKAGES	########################
library(dada2)
library(ShortRead)
library(tidyverse)
library(biomformat)

############	READ THE FILES AND THE ENV VARIABLES	############

# Get environmet variable
# PHRED SCORE
PHRED<-Sys.getenv("PHRED")
# PATH TO SCRIPTS
DIR<-Sys.getenv("REF_DATA_PATH")
# Number of cores to use
NCORES<-Sys.getenv("NSLOTS")
NCORES<-as.integer(NCORES)
print(paste("Number of cores used",NCORES))
# Set the working directory
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd(getSrcDirectory()[1])

# Read list of final names
#final_names<-read_tsv("Final_name.txt", col_names = F)
#metatable<-read_tsv("Metatable.txt", col_names = T)
final_names<-read_tsv("Final_name.txt", col_names = F)
metatable<-read_tsv(Sys.getenv("METATABLE"), col_names = T)
# Select only the sample ID and run columns
metatable<-select(metatable, matches("#SampleID|run"))
# Left join final_name and metatable by sample ID
# Change "_" to "-" to match the Metatable 
final_names$X2<-gsub("_","-",final_names$X2)
final_names<-left_join(final_names,metatable, by = c("X2" = "#SampleID"))
colnames(final_names)[ncol(final_names)]<-"X3"


#############################################################################################################################
################	RUN DADA2 IN PAIRED-END MODE	#####################################################################

# Sickle folder name
sickle_folder<-paste("Q",PHRED,sep = "")
# Get paths to fastq
fastq_path<-list()
x<-1
for(directory in unique(final_names$X2)){
  fastq_path[[x]]<-file.path(paste(directory,"Cutadapt_Sickle",sickle_folder, sep="/"))
  x<-x+1
}
fastq_path<-unlist(fastq_path)
# Get fastq names 
forward_names<-list()
reverse_names<-list()
x<-1
for(directory in fastq_path){
  forward_names[[x]]<-list.files(directory, pattern="_R1")
  reverse_names[[x]]<-list.files(directory, pattern="_R2")
  x<-x+1
}
forward_names<-unlist(forward_names)
reverse_names<-unlist(reverse_names)
# Name of filterd fastq files
x<-1
filtFs<-list()
filtRs<-list()
for(directory in unique(final_names$X2)){
  filtFs[[x]] <- paste0(unique(final_names$X2)[x], "_F_DADA2filt.fastq.gz")
  filtRs[[x]] <- paste0(unique(final_names$X2)[x], "_R_DADA2filt.fastq.gz")
  x<-x+1
}
filtFs<-unlist(filtFs)
filtRs<-unlist(filtRs)

print("Trimming reads in paired-end mode")
# Run the trimming in parallel mode
filterAndTrim(fwd=file.path(fastq_path, forward_names), filt=file.path(fastq_path, filtFs),
              rev=file.path(fastq_path, reverse_names), filt.rev=file.path(fastq_path, filtRs),
              truncLen=c(0,0), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=NCORES, matchIDs=TRUE)


# Sample inference and merger of paired-end reads
mergers <- vector("list", length(unique(final_names$X2)))
names(mergers) <- unique(final_names$X2)
names(filtFs) <- unique(final_names$X2)
names(filtRs) <- unique(final_names$X2)
names(fastq_path)<-unique(final_names$X2)

# The following learns the error rates only once for each run
# It requires a third column in the final_names specifying the run for each sample
for(runs in unique(final_names$X3)){
	cat("Processing:", "Analysing Samples in Run:",runs, "\n")
	# Filter to only the samples of the current run
	final_names_sub<-final_names %>% filter(X3 == runs)
	# Filter the forward, reverse and paths for samples sequenced in the same run
	filtFs_run <- filtFs[unique(final_names_sub$X2)]
	filtRs_run <- filtRs[unique(final_names_sub$X2)]
	fastq_path_run <- fastq_path[unique(final_names_sub$X2)]
	
	cat("Processing:", "Learning errors for the current run","\n")
	# Learn errors for the current run
	cat("Forward","\n")
	errF <- learnErrors(paste(fastq_path_run,filtFs_run, sep="/"), nbases=1e8, multithread=NCORES)
	cat("Reverse","\n")
	errR <- learnErrors(paste(fastq_path_run,filtRs_run, sep="/"), nbases=1e8, multithread=NCORES)
	# Plot error rates by Quality Score (Phred)
	p1<-plotErrors(errF, nominalQ=TRUE)
	ggsave(paste("Estimated_error_rates_ForwardReads_",runs,".pdf", sep=""), plot=p1, device="pdf", path="DADA2_OTU_tables/Error_rate_plots")
	# Plot error rates by Quality Score (Phred)
	p2<-plotErrors(errR, nominalQ=TRUE)
	ggsave(paste("Estimated_error_rates_ReverseReads_",runs,".pdf", sep=""), plot=p2, device="pdf", path="DADA2_OTU_tables/Error_rate_plots")
	
	for(nam in unique(final_names_sub$X2)){
		cat("Processing:", nam, "\n")
		#Dereplicate the filtered forward fastq reads
		cat("Processing:", nam,"Dereplicate forward", "\n")
		derepF <- derepFastq(paste(fastq_path[[nam]],filtFs[[nam]], sep="/"))
		# Run dada on forward
		cat("Processing:", nam,"Run DADA2 on dereplicated forward reads", "\n")
		ddF <- dada(derepF, err=errF, multithread=NCORES)
  
		#Dereplicate the filtered reverese fastq reads
		cat("Processing:", nam,"Dereplicate reverse", "\n")
		derepR <- derepFastq(paste(fastq_path[[nam]],filtRs[[nam]], sep="/"))  
		# Run dada on reverse
		cat("Processing:", nam,"Run DADA2 on dereplicated reverse reads", "\n")
		ddR <- dada(derepR, err=errR, multithread=NCORES)
  
		# Merge sequences
		cat("Processing:", nam,"Merging", "\n")
		merger <- mergePairs(ddF, derepF, ddR, derepR)
		mergers[[nam]] <- merger
	}
}

####################################################################################
#############################	 WRITE THE RESULTS	############################

print("Make table")
#Construct sequence table
seqtab <- makeSequenceTable(mergers)
print("Remove chimera")
#Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)


print("Assign taxonomy")
#Assign taxonomy:
taxa <- assignTaxonomy(seqtab.nochim, paste(DIR,"SILVA/DADA2_SILVA123/silva_nr_v123_train_set.fa.gz", sep="/"))
taxa.plus <- addSpecies(taxa, paste(DIR,"SILVA/DADA2_SILVA123/silva_species_assignment_v123.fa.gz", sep="/"), verbose=TRUE, tryRC=T, allowMultiple=TRUE)
head(unname(taxa.plus))

colnames(taxa.plus) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")

print("Save sequence file: DADA2_OTU_tables/seqs.fa")
#Save sequence file:
sink("DADA2_OTU_tables/seqs.fa");cat(paste(">","SEQ_",seq(1:dim(seqtab.nochim)[2]),"\n",paste(colnames(seqtab.nochim),"\n",sep=""),sep=""),sep="");sink()

seqtab.final<-seqtab.nochim
colnames(seqtab.final)<-paste("SEQ_",seq(1:dim(seqtab.nochim)[2]),sep="")

seqtab.final<-t(seqtab.final)
seqtab.final<-as.data.frame(seqtab.final)
seqtab.final$SEQ<-rownames(seqtab.final)
seqtab.final<-seqtab.final[,c(ncol(seqtab.final),1:ncol(seqtab.final)-1)]

print("Save sequence table: DADA2_OTU_tables/seq_table.txt")
#Save the sequence table:
write.table(seqtab.final,"DADA2_OTU_tables/seq_table.txt",quote=F, sep="\t", row.names=F)

taxa.plus.final<-as.data.frame(apply(unname(taxa.plus), 1, paste, collapse=";"))
taxa.plus.final$OTUID<-paste("SEQ_",seq(1:dim(seqtab.nochim)[2]),sep="")
taxa.plus.final<-taxa.plus.final[,c(2,1)]
colnames(taxa.plus.final)<-c("OTUID","taxa")
write.table(taxa.plus.final, "DADA2_OTU_tables/seq_Taxonomy.txt",sep="\t",na="", row.names=F, quote=F)


