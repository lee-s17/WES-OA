R version 4.3.1 (2023-06-16) -- "Beagle Scouts"

setwd("~/Documents/workingdir/WES/jehai_temuan/postVC/e110post")
#input matchvep file
library(dplyr)
tab <- read.table("matchvep-cut-ORG-annpass-OA_snp.vcf", header = TRUE)
snp <- tab %>% distinct() #remove duplicate lines

#####Variant Consequences
#get variants with only one result-unique, or more than one conseq
unqvep <- snp %>%
  group_by(Uploaded_variation) %>%
  filter(n()==1)
unqvep_table <- data.frame(table(number=unqvep$Consequence))
#check variants with more than one conseq
dupvep <- snp %>% 
  group_by(Uploaded_variation) %>% 
  filter(n()>1)

all_vep <- rbind(unqvep, dupvep)
write.table(all_vep, "conseq-matchvep-cut-ORG-annpass-OA_snp", quote = FALSE, row.names = FALSE, sep = "\t")



#####Ancestral allele
ancestral <- select(snp, 1,2,37) %>% distinct()
all_anc_table <- data.frame(table(number=ancestral$AA)) #get count of variants with no AA

bcflist <- read.table("bcfann-ORG-annpass-OA_snp.txt", header = TRUE)
bcflist$chrompos <- paste0(bcflist$CHROM,":",bcflist$POS)
