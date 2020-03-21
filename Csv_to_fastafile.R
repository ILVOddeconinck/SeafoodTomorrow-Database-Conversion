#######################
#Required packeges

library("seqRFLP")
library("Biostrings")

#Note, for this you need to have Ncbi blast+ installed.
#See: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download For more information

#Remember to use setwd() in case the files you need are not within the project folder
#You can also use ctrl+shift+H to select the folder you work in 

#Read in the excel file with SFT sequences
data = read.csv(file = "Reference species.csv", header = TRUE, sep = ";", quote = "\"", dec = ".", fill = TRUE, comment.char = "",fileEncoding = "UTF-16LE")

#The dataset needs to be cleaned up. In the following steps we will slowly transform the set into a fasta file.

#Create a data frame that only contains teh essentials.
data = data.frame(data$Reference.species.Scientific.name, data$Specimen.Code, data$DNA.Sequence.GeneID, data$DNA.Sequence.DNA.sequence)
colnames(data) = c("Scientific_name", "Specimen_code", "Gene", "Sequence")

#Transform the sequences to readable sequences by removing characters that aren't nucleotides
data$Sequence = gsub(" ", "", data$Sequence)
data$Sequence = gsub("\\|", "", data$Sequence)

#Replace the space in Scientific name with an underscore. Spaces are problematic later on in the script
data$Scientific_name = gsub(" ", "_", data$Scientific_name)

#Create a sequence ID for our fasta file. This contains: Scientific name, specimen code and Gene
SequenceID = paste(data$Scientific_name,"_",data$Specimen_code,"_",data$Gene, sep = "")

#Create a database with the sequence ID and the Sequences.
Database = data.frame(SequenceID, data$Sequence)
colnames(Database) = c("SequenceID", "Sequence")

#Remove empty rows
Database = Database[!(Database$Sequence == ""), ]

#Turn this database into a fasta format
df.fasta = dataframe2fas(Database, file = "Database.txt")
DatabaseFasta = readDNAStringSet("Database.txt")

#Write the FASTA file into the map. Note that you can put this anywhere you want.
writeXStringSet(DatabaseFasta,filepath = "SeafoodtomorrowDB.fasta")
