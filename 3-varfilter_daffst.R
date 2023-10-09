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
SNPVEP_BCF <- merge(x=snpvep, y=bcfsnp, by.x ="Location", by.y = "chrompos", all.y =  TRUE) %>%
  dplyr::rename("AF_1KGP"=AF.x,"AA_VEP.e110"=AA,"AF_WES"=AF.y)

wes_hgdp <- read.table("/home/z6/Documents/workingdir/WES/postVC/e110post/hgdp/wes-hgdporg", header=TRUE)
missing_hgdp <- read.table("/home/z6/Documents/workingdir/WES/postVC/e110post/hgdp/bcfannfreq-missinghgdp", header=TRUE)
wes_all <- rbind(wes_hgdp,missing_hgdp)
wes_all$chrompos <- paste0(wes_all$CHROM,":",wes_all$POS)
wes_all <- select(wes_all, c(3:7,11,14,16)) %>% rename("REF_HGDP"=REF,"ALT_HGDP"=ALT,"AN_HGDP"=AN, "ID_HGDP"=ID)
i <- c(6, 7)
wes_all[ , i] <- apply(wes_all[ , i], 2, function(x) as.numeric(as.character(x)))
sapply(wes_all, class)
SNPVEP_BCF_HGDP <- merge(x=SNPVEP_BCF, y=wes_all, by.x ="Location", by.y = "chrompos", all.x =  TRUE)
SNPVEP_BCF_HGDP<-  select(SNPVEP_BCF_HGDP, !c(12,15,20))



#####dDAF and Fst
ancestralSNPVEP_BCF_HGDP <- select(SNPVEP_BCF_HGDP, 1,2,34,40:42,46,49,51:57)  %>% distinct() 
ancestralSNPVEP_BCF_HGDP <- ancestralSNPVEP_BCF_HGDP %>%
  mutate(DAF_JH = ifelse(REF == toupper(AA_VEP.e110), AF_JH, ifelse(ALT == toupper(AA_VEP.e110), 1 - AF_JH,
                                                                    ifelse(REF == AA_chimp, AF_JH, ifelse(ALT == AA_chimp, 1 - AF_JH, NA))))) %>%
  mutate(DAF_TM = ifelse(REF == toupper(AA_VEP.e110), AF_TM, ifelse(ALT == toupper(AA_VEP.e110), 1 - AF_TM,
                                                                    ifelse(REF == AA_chimp, AF_TM, ifelse(ALT == AA_chimp, 1 - AF_TM, NA))))) 
ancestralSNPVEP_BCF_HGDP <- ancestralSNPVEP_BCF_HGDP %>%
  mutate(DAF_AFR = ifelse(REF_HGDP == toupper(AA_VEP.e110), AF_AFR , ifelse(ALT_HGDP == toupper(AA_VEP.e110), 1 - AF_AFR,
                                                                    ifelse(REF_HGDP == AA_chimp, AF_AFR, ifelse(ALT_HGDP == AA_chimp, 1 - AF_AFR, NA))))) %>%
  mutate(DAF_NONAFR = ifelse(REF_HGDP == toupper(AA_VEP.e110), AF_NONAFR, ifelse(ALT_HGDP == toupper(AA_VEP.e110), 1 - AF_NONAFR,
                                                                    ifelse(REF_HGDP == AA_chimp, AF_NONAFR, ifelse(ALT_HGDP == AA_chimp, 1 - AF_NONAFR, NA))))) 

ancestralSNPVEP_BCF_HGDP$dDAF_WES<- abs(ancestralSNPVEP_BCF_HGDP$DAF_JH- ancestralSNPVEP_BCF_HGDP$DAF_TM)
ancestralSNPVEP_BCF_HGDP$dDAF_WES<- as.numeric(as.character(format(round(ancestralSNPVEP_BCF_HGDP$dDAF_WES, 4), nsmall = 4)))
ancestralSNPVEP_BCF_HGDP$dDAF_HGDP<- abs(ancestralSNPVEP_BCF_HGDP$DAF_AFR- ancestralSNPVEP_BCF_HGDP$DAF_NONAFR)
ancestralSNPVEP_BCF_HGDP$dDAF_HGDP<- as.numeric(as.character(format(round(ancestralSNPVEP_BCF_HGDP$dDAF_HGDP, 4), nsmall = 4)))
                      
fst <- read.table("fstWES-ORG-annpass-genohweking_biallsnp.JH.TM.fst.var",header = TRUE)
fstx <- read.table("fstWES-ORG-annpass-genohweking_biallsnp.x.JH.TM.fst.var",header = TRUE)
fstall <- rbind(fst,fstx)
fstall$chrompos <- paste0("chr",fstall$ID)

SNPVEP_BCF_HGDP_ANC <- merge(x=SNPVEP_BCF_HGDP, y=select(ancestralSNPVEP_BCF_HGDP, 1,16:21), by.x ="Location", by.y = "Location", all.x = TRUE)
SNPVEP_BCF_HGDP_ANC_FST <- merge(x=SNPVEP_BCF_HGDP_ANC, y=select(fstall, 4:6), by.x ="Location", by.y = "chrompos", all.x = TRUE)
write.table(SNPVEP_BCF_HGDP_ANC_FST, "rfile/SNPVEP_BCF_HGDP_ANC_FST.txt", quote = FALSE, row.names = FALSE, sep = "\t")


#quantile and outliers
ddaf <- filter(ancestralSNPVEP_BCF_HGDP,ancestralSNPVEP_BCF_HGDP$AN >= 90 & !is.na(ancestralSNPVEP_BCF_HGDP$AF_JH )) #genotype rate
quantile(ddaf$dDAF_WES, c(.99), na.rm = TRUE)
upperdaf <- filter(ddaf,ddaf$dDAF_WES> 0.37438) 
breaks <- c(-0.1, 0.1, 0.2 ,0.3, 0.4,0.5)
pivotdaf <- cut(ddaf$dDAF_WES, breaks= breaks, labels = NULL)
daf <- data.frame(Range=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","NA"),
                  Value="dDAF", SNPs=c(summary(pivotdaf)))

hfst <- select(filter(SNPVEP_BCF_HGDP_ANC_FST,SNPVEP_BCF_HGDP_ANC_FST$AN >= 90) ,c(1,2,65))%>% distinct() 
quantile(hfst$HUDSON_FST, c(.99), na.rm = TRUE)
upperfst <- filter(hfst,hfst$HUDSON_FST> 0.2570786) 
pivotfst <- cut(hfst$HUDSON_FST, breaks= breaks, labels = NULL)
fst <- data.frame(Range=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5","NA"),
                  Value="FST", SNPs=c(summary(pivotfst)))          
#plot
library(ggplot2)
ggplot(rbind(daf,fst), aes(x = factor(Range), y = SNPs, fill = Value, colour = Value)) +
  geom_bar(stat = "identity", position = "dodge")+theme_bw(base_size = 30)+labs(y= "No. of SNPs", x = "Bin")+
  geom_text(aes(label = SNPs), vjust = -0.5, colour = "black",size=8,position = position_dodge(.9))
