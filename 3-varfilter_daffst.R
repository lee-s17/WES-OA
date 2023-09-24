R version 4.3.1 (2023-06-16) -- "Beagle Scouts"

setwd("~/Documents/workingdir/WES/jehai_temuan/postVC/e110post")
#input matchvep file
library(dplyr)
tab <- read.table("matchvep-cut-ORG-annpass-OA_snp.vcf", header = TRUE)
snp <- tab %>% distinct() #remove duplicate lines



#####Variant Consequences
#get variants with only one result-unique, or more than one conseq
unqvep <- snp %>% group_by(Uploaded_variation) %>%
  filter(n()==1)
unqvep_table <- data.frame(table(number=unqvep$Consequence))
#check variants with more than one conseq
dupvep <- snp %>% group_by(Uploaded_variation) %>% 
  filter(n()>1)
all_vep <- rbind(unqvep, dupvep)
write.table(all_vep, "conseq-matchvep-cut-ORG-annpass-OA_snp", quote = FALSE, row.names = FALSE, sep = "\t")



#####Compile masterfile: bcftools, vep, hgdp
snpvep <- filter(snp, CANONICAL =="YES")
bcfsnp <- read.table("bcfann-ORG-annpass-OA_snp.txt", header = TRUE, stringsAsFactors=FALSE) 
sapply(bcfsnp, class)
i <- c(10, 13)
bcfsnp[ , i] <- apply(bcfsnp[ , i], 2, function(x) as.numeric(as.character(x)))

bcfsnp$chrompos <- paste0(bcfsnp$CHROM,":",bcfsnp$POS)
SNPBCF_VEP <- merge(x=snpvep, y=bcfsnp, by.x ="Location", by.y = "chrompos", all.y =  TRUE) %>%
  dplyr::rename("AF_1KGP"=AF.x,"AA_VEP.e110"=AA,"AF_WES"=AF.y)

hgdp <- read.table("/home/z6/Documents/workingdir/hgdpdec20/hgdp/final/daf/dafedit_afr.non.txt", header=TRUE)
hgdpfst <- read.table("/home/z6/Documents/workingdir/hgdpdec20/hgdp/final/fst/fst-afrnon-orgset-genohwekingbiallsnpsauto.AFR.NON.fst.var", header=TRUE)
hgdpall <- merge(x=hgdp, y=hgdpfst, by.x = "CHROMPOS", by.y = "ID", all.x =  TRUE)  
hgdpall  <-  select(hgdpall, !c(2,10,11,15,16))%>% dplyr::rename("CHROM"=CHROM.x, "POS"=POS.x, "REF_HGDP"=REF.x,"ALT_HGDP"=ALT.x,"AA_HGDP"=AA, "AA_HGDP"=AA)
hgdpall$chrompos <- paste0("chr",hgdpall$CHROMPOS)

SNPBCF_VEP_HGDP <- merge(x=SNPBCF_VEP, y=hgdpall, by.x ="Location", by.y = "chrompos", all.x =  TRUE)
SNPBCF_VEP_HGDP<-  select(SNPBCF_VEP_HGDP, !c(55,56,57))



#####dDAF and Fst
ancestralSNPBCF_VEP <- select(SNPBCF_VEP_HGDP, 1,37,40,41,42,43,44,49,52,55,56,57,59:63)  %>% distinct() 
ancestralSNPBCF_VEP <- ancestralSNPBCF_VEP %>%
  mutate(DAF_JH = ifelse(REF == toupper(AA_VEP.e110), AF_JH, 1 - AF_JH )) %>%
  mutate(DAF_TM = ifelse(REF == toupper(AA_VEP.e110), AF_TM, 1 - AF_TM ))
ancestralSNPBCF_VEP$dDAF_WES<- abs(ancestralSNPBCF_VEP$DAF_JH- ancestralSNPBCF_VEP$DAF_TM)
ancestralSNPBCF_VEP$dDAF_WES<- as.numeric(as.character(format(round(ancestralSNPBCF_VEP$dDAF_WES, 4), nsmall = 4)))

fst <- read.table("fstWES-ORG-annpass-genohweking_biallsnp.JH.TM.fst.var",header = TRUE)
fstx <- read.table("fstWES-ORG-annpass-genohweking_biallsnp.x.JH.TM.fst.var",header = TRUE)
fstall <- rbind(fst,fstx)

SNPBCF_VEP_HGDP_ANC <- merge(x=SNPBCF_VEP_HGDP, y=select(ancestralSNPBCF_VEP, 1,18,19,20), by.x ="Location", by.y = "Location", all.x = TRUE)
SNPBCF_VEP_HGDP_ANC_FST <- merge(x=SNPBCF_VEP_HGDP_ANC, y=select(fstall, 3,4,5), by.x ="CHROMPOS", by.y = "ID", all.x = TRUE)

quantile(SNPBCF_VEP_HGDP_ANC_FST$dDAF_WES, c(.99), na.rm = TRUE)
quantile(SNPBCF_VEP_HGDP_ANC_FST$HUDSON_FST, c(.99), na.rm = TRUE)
breaks <- c(-0.1, 0.1, 0.2 ,0.3, 0.4,0.5)
# bucketing values into bins
pivot <- cut(SNPBCF_VEP_HGDP_ANC_FST$dDAF_WES, breaks= breaks, labels = NULL)
data.frame(dDAF=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","NA"),
           dDAF_WES=c(summary(pivot)))
upperdaf <- filter(SNPBCF_VEP_HGDP_ANC_FST,SNPBCF_VEP_HGDP_ANC_FST$dDAF_WES> 0.3734) 
upperfst <- filter(SNPBCF_VEP_HGDP_ANC_FST,SNPBCF_VEP_HGDP_ANC_FST$HUDSON_FST> 0.2567) 
