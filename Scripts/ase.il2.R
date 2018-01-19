library(data.table)
library(XLConnect)
library(biomaRt)
library(lme4)

source('/home/ev250/ase-il2ra/Functions/various.R')
source('/home/ev250/Bayesian_inf/trecase/Functions/various.R')  # for hap description

## read Dan Rainbow's data

wb <- loadWorkbook("/mrc-bsu/scratch/ev250/ase-il2ra/objects/CM_ASE_genotypes.xlsx")

## extract relevant info

DT <- data.table(readWorksheet(wb,sheet="Tag SNPs",startRow=5,endRow=69,startCol=1,endCol=53)) # genotype tag SNPs
DT2 <- data.table(readWorksheet(wb,sheet="Tag SNPs",startRow=5,endRow=69,startCol=1,endCol=55)) # genotype tag SNPs and ASE

haps <- unique(data.table(t(readWorksheet(wb,sheet="Tag SNPs",startRow=4,endRow=5,startCol=2,endCol=53,header=F)))) # indicates to which group each SNP belongs to.
names(haps) <- c("hap","SNP")
haps[grep("rs",hap),hap:="Sus"]
haps[,protective:=c("A","A","T","G","T","G","C","A","C","G","T","A","C","G","T","T","A","G","T","T","T","G","A","T","C","A")]
haps[, susceptible:=c("G","G","A","T","C","A","T","G","T","A","C","T","T","A","C","C","G","T","C","C","C","A","G","C","T","G")]

ASE <- data.table(readWorksheet(wb,sheet="ASE data",startRow=1,endRow=65,startCol=1,endCol=3)) # ASE data is not ratio G:A but percetange of G allele.
## sort ASE in the same order as DT$ID

#ASE[match(DT$ID,Donor),][!is.na(Donor)]

## work out which haplotype each individual has: follow Daniel's rules:

## A = rs61839660 is T (there are some donors that have recombined the A haplotype, but as we are pretty sure rs61839660 is the causal variant for the A haplotype, then anyone carrying the T protective allele is counted as A haplotype).
## A-660 = rs61839660 is C and rs12722495 is C
## C = rs11594656 is A and rs11597367 is G
## D = rs56382813 is T and rs41295055 is T and rs3118475 is C
## E = rs6602437
## F = is very rare and only one or two people carry it and i have been using rs41295121 as T
## Susceptible is classed as not carrying any of the A, C, D or E alleles.

haps.D <- haps[SNP %in% c("rs61839660", "rs12722495",  "rs11594656","rs11597367", "rs56382813","rs41295055","rs3118475", "rs6602437", "rs41295121"), ]
                          
DT <- m.hap(haps.D,DT)

## added haplotype cols to DT
hap.cols <- unlist(lapply(unique(haps$hap), function(i) grep(paste0("^",i), x=names(DT), value=T)))

## select het individuals for hap A

hets <- DT[A.1!=A.2 & A.2=="A.s" ,c("ID",hap.cols), with=F]

## recode haplopytes: A, A-, C, D, E with 0 (homo susceptible) or 1 (het)

hets.r <- r.hap(hets,groups=c("A","A-","C","D","E"), haps.D) # recode A hets selecting the same hap groups as in Dan's report

## add ASE

hets.ASE <- merge(hets.r,ASE, by.x="ID", by.y="Donor", all.x=T)

## remove NA

hets.ASE <- hets.ASE[complete.cases(hets.ASE),]

## add col percetage A_allele (ASE) 100-%G

hets.ASE[, p.A:=100-rs12244380.G.A]

#############################################################################################################
#################################### Phasing ################################################################

## prepare BED file for lifting B36 coordinates into B37 for phasing

coord <- data.table(t(readWorksheet(wb,sheet="All ichip SNPs",startRow=1,endRow=4,startCol=3, endCol=1023, header=F)))
names(coord) <- unlist(coord[1,])
coord <- coord[-1,]
coord[,`Build 36 position`:=as.numeric(`Build 36 position`)]

##BED file uses pos-1 (first base of chr=0)

bed <- coord[,.(Chr,`Build 36 position`, SNP)]
bed[,Chr:=paste0("chr",Chr)][,`Build 36 position`:=`Build 36 position`-1][is.na(SNP),SNP:=paste0(Chr,":",`Build 36 position`)]
bed[,pos.end:=`Build 36 position`]

## save
write.table(unique(bed[,c(1:2,4),with=F]), file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/ichipB36.bed",row.names=F,col.names=F,quote=F)

## Convert to B37 using CrossMap run in bash

module load python/2.7.5
source /scratch/ev250/bin/PYTHON/bin/activate
/scratch/ev250/bin/PYTHON/bin/CrossMap.py bed \
/mrc-bsu/scratch/ev250/reference_genome/36to37/NCBI36_to_GRCh37.chain.gz \
/mrc-bsu/scratch/ev250/ase-il2ra/objects/ichipB36.bed \
/mrc-bsu/scratch/ev250/ase-il2ra/objects/ichipB37.bed

## Make PED/MAP files for phasing

geno <-  data.table(readWorksheet(wb,sheet="All ichip SNPs",startRow=6,endRow=57,startCol=1, endCol=1023, header=F))
## remove col 3, not clear what it is from Dan's file but is 1 for all samples
geno[,Col3:=NULL]

coord[is.na(SNP),SNP:=paste0(Chr,":",`Build 36 position`)] # replace NA with chr:pos b36
names(geno) <- c("ID","sex",paste0(coord$SNP, c(".1",".2")))

##format geno as PED
ped <- geno[,FamID:=ID][,c("FID","MID","Ph"):= double(nrow(geno))]
##order as ped
setcolorder(ped, c("FamID","ID","FID","MID","sex","Ph",paste0(coord$SNP, c(".1",".2")) ))

## remove snps with only missing values (all 0 entries)
s <- seq(7,ncol(ped),2)

ss <- sapply(s, function(i) sum(unique(unlist(ped[,i:(i+1),with=F]))==0)==1 & length(unique(unlist(ped[,i:(i+1),with=F]))==0)==1) # length==1 excludes some missing values TRUE,FALSE vector
mis.snps <- names(ped)[c(s[which(ss==TRUE)], s[which(ss==TRUE)]+1)] ## to account for .1 and .2
ped[,(mis.snps):=NULL]

## identify snps with I/D and remove, Shape it doesnt take I/D for phasing
m.v <-apply(geno,2,function(i) which(sum(i=="I" | i=="D")>0))

m.snps <- unique(names(m.v[unlist(lapply(m.v, function(i) length(i)>0))])) # get snps with insertions/del
ped[,(m.snps):=NULL]

##  add b37 coordinates to coord
b37 <- fread("/mrc-bsu/scratch/ev250/ase-il2ra/objects/ichipB37.bed")
u.coord <-cbind(unique(coord), b37=b37$V2+1) ## add 1 because it is bed format

##map file : chr, SNPID, genetic pos (0) and b37 pos

## remove m.snps &  mis.snps from u.coord

u.noindels <- u.coord[!SNP %in% unique(gsub("\\..*","",c(m.snps,mis.snps))),]
u.noindels[,gendis:=0]


## use biomart to check all sites, also need to update snp names as in reference panel
bmart <- useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp", host="grch37.ensembl.org", path="/biomart/martservice")
snpcheck = unique(rbindlist(lapply(1:nrow(u.noindels), function(i) data.table(getBM(attributes = c('refsnp_id','allele','chrom_start'), filters = c('chr_name','start','end'), values=list(u.noindels[i,Chr], u.noindels[i,b37], u.noindels[i,b37] +1) ,mart=bmart)))))

## remove indels from snpcheck

idls <- which(sapply(strsplit(gsub("/","",snpcheck$allele),""), function(i) "-" %in% i))
snpcheck <- snpcheck[-idls,]

## merge snpcheck with u.noindels

u.noind <- merge(u.noindels,snpcheck,by.x="b37",by.y="chrom_start", all.x=T)

## get alleles coded by Dan excluding missing values (0)
all.dan <- sapply(seq(7,ncol(ped),2),function(i) {
    u=unique(unlist(ped[,i:(i+1)]),sep="/")
    u=u[which(u!="0")] ## exclude missing value "0"
    return(u)
    })      
names(all.dan) <- unique(gsub("\\..*","", names(ped)[7:ncol(ped)]))

## match all.dan snps from u.noind : 
s <- c()
for(i in seq_along(all.dan)){
    tmp <- strsplit(u.noind[SNP %in% names(all.dan)[[i]], allele], "/")
    
    s <- c(s, paste0(unlist(tmp[which(sapply(tmp, function(j) sum(all.dan[[i]] %in% j))==length(all.dan[[i]]))]), collapse="/"))
}


## save map with SNP names from biomart to match reference panel

map <- u.noind[allele %in% s,]

write.table(map[,.(Chr,refsnp_id,gendis,b37)], file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/snpsB37no_I_D.map",row.names=F,col.names=F,quote=F)

## save it also with all columns
write.table(map, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/snpsB37no_I_D_map.txt",row.names=F)


## replace SNP names in ped with updated names refsnp_id

for(i in 1:nrow(map)){
    g <- grep(map$SNP[i],names(ped))
    names(ped)[g] <- paste0(map$refsnp_id[i], c(".1",".2"))
}

## same samples only have genotype info for tags, add info for these samples into ped file
geno.tag <- data.table(readWorksheet(wb,sheet="Tag SNPs",startRow=5,endRow=69,startCol=1,endCol=55))
## recode snp names in geno.tag
names(geno.tag) <- c("ID", paste0(gsub("_[1,2]","",gsub("\\.1","",names(geno.tag)[2:ncol(geno.tag)])), c(".1",".2")))

## identify samples with missing data
samp.m <- geno.tag[!ID %in% ped$ID,]

## add samp.m to ped file, coding missing data with 0

## prepare ped file from samp.m, get missing snps

m.s <- names(ped)[!names(ped) %in% names(geno.tag)]
ped.m <- data.table(matrix(0,nrow=nrow(samp.m),ncol=ncol(ped)))
names(ped.m) <- names(ped)
ped.m[, names(samp.m):=lapply(1:ncol(samp.m), function(i) unlist(samp.m[,i,with=F]))] 
ped.m[,FamID:=ID]

## append ped.m to ped

ped <- rbind(ped,ped.m)

## save ped file

write.table(ped, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/genoB37no_I_D.ped",row.names=F,col.names=F,quote=F)

## also as txt with col names
write.table(ped, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/genoB37no_I_D.txt", row.names=F)


## run shapeit using shapeit_no_indels.sh
## read shapeit -check output file: SNPs with aligments problems
## aligment problems: homozygous individuals from ped file give a wrong alt allele (ref allele) and some snps are missing in ref panel

##check reference panel

rp <- fread("zcat '/scratch/wallace/1000GP_Phase3/1000GP_Phase3_chr10.legend.gz'")

rp <- rp[position %in% map$b37,]

alig.prob <- fread("/mrc-bsu/scratch/ev250/ase-il2ra/shapeit/chr10.aligments.snp.strand")

missing.in.fp <- alig.prob[V1=="Missing",] ## snps truly missing in reference panel

## select snps only homo in my samples excluding indels

homo <- alig.prob[V1=="Strand" & nchar(V9)==1 & nchar(V10)==1,]

########################################## Phasing second strategy ##################################

### try phasing with an extra faked sample het for the homo snps so I can use that info when phasing and see if that makes any difference

l1 <- lapply(c(rep("fake",2), rep(0,(ncol(ped)-2))), function(i) i) ## make list to append to ped
ped.f <- rbindlist(list(ped, l1))

## add geno info for homo snps to fake sample
h.snps <- sapply(homo[,V4], function(i) paste0(i, c(".1",".2")))

v <- as.vector(t(as.matrix(homo[,V9:V10])))

for(i in seq_along(h.snps)) {set(ped.f, nrow(ped.f) , h.snps[i] , v[i]) } ## replace values in row

## save ped file

write.table(ped.f, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/genoB37no_I_D_hetsample.ped",row.names=F,col.names=F,quote=F)

## also as txt with col names
write.table(ped.f, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/genoB37no_I_D_hetsample.txt", row.names=F)


## convert ped.f and map into bed/bim/fam for phasing so ref and alt allele are defined, then remove fake sample when phasing

##in bash, module load plink
## cd /mrc-bsu/scratch/ev250/ase-il2ra/objects/

plink --ped genoB37no_I_D_hetsample.ped --map snpsB37no_I_D.map --make-bed --out genoB37no_I_D_hetsample --missing-genotype 0


########################### Phasing strategy 3 (didnt run it) ####################################
## create GEN/sample files using gtool

/mrc-bsu/scratch/ev250/bin/gtool -P --ped genoB37no_I_D_hetsample.ped --map snpsB37no_I_D.map --og genoB37no_I_D_hetsample.gen --os genoB37no_I_D_hetsample.sample

## files are OK but as with bed/bim/fam the selection of which allele is rf and which one is alt is unpredictable, and different from bim/bed/fam.

######################################################################################
################################# Post-phasing ######################################

samples <- fread("cat /mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.phased.with.ref.sample | cut -d' ' -f1 | tail -n+3", header=F) ## remove first 2 rows

haps.sam <- fread("/mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.phased.with.ref.haps", header=F)

## select Dan snps: c(haps.D$SNP, "rs12244380") from haps

snps.d <- c(haps.D[ ,SNP], "rs12244380")

haps.samp <- haps.sam[V2 %in% snps.d,]

## add phasing info to samples and recode haps using Dan's rules for haps A,A-660,C,D and E: 0 susceptible and 1 protective, 9 codes for missing value

haps.D <- rbind(haps.D, list("ASE", "rs12244380", "A","G"))

samples <- phased.r(samples,haps.samp, haps.D)

## select samples with ASE info
## recode "ASE" Donor to uppercase and the same for "samples" V1
samples[,V1:=toupper(V1)]          
ASE[,Donor:=toupper(Donor)]

sam.AS <- merge(ASE,samples, by.x="Donor", by.y="V1", all=T)



############################################ strategy 2 ############################################

### A using ped.f and map for phasing #########################

samples2 <- fread("cat /mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.het.phased.with.ref.sample | cut -d' ' -f1 | tail -n+3", header=F) ## remove first 2 rows

haps.samp2 <- fread("/mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.het.phased.with.ref.haps", header=F)
haps.samp2 <- haps.samp2[V2 %in% snps.d,]
samples2 <- phased.r(samples2,haps.samp2, haps.D)

## remove fake sample
samples2 <- samples2[V1!="fake",]

## add ASE, recode "ASE" Donor to uppercase and the same for "samples" V1
samples2[,V1:=toupper(V1)]          
sam.AS2 <- merge(ASE,samples2, by.x="Donor", by.y="V1", all=T)


##### B using bed bim fam for phasing and excluding "fake" sample  when running shapit #########################
samples3 <- fread("cat /mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.bbfam.phased.with.ref.sample | cut -d' ' -f1 | tail -n+3", header=F) ## remove first 2 rows

haps.samp3 <- fread("/mrc-bsu/scratch/ev250/ase-il2ra/shapeit/il2ra.noindels.bbfam.phased.with.ref.haps", header=F)
haps.samp3 <- haps.samp3[V2 %in% snps.d,]


samples3 <- phased.r(samples3,haps.samp3, haps.D)

## add ASE, recode "ASE" Donor to uppercase and the same for "samples" V1
samples3[,V1:=toupper(V1)]          
sam.AS3 <- merge(ASE,samples3, by.x="Donor", by.y="V1", all=TRUE)


###############################################################################################
## QC of phased samples: select snps.d from ped and check if genotype is compatible with haps.samp with ped

## format haps.samp with first col ID and then snp1.1, snp1,2, etc as in ped
h.ped <- haps.ped(samples,haps.samp)

## recode ped as haps.samp (0,1)
ped.d <- ped[,names(ped)[which(names(ped) %in% names(h.ped)[2:ncol(h.ped)])],with=F]
ped.d <- r.ped(ped.d,haps.samp)

## check consistency, for each snp proportion of concordance between genotypes

cons <- h.p.comp(ped.d,h.ped,snps.d)


####################################################################################
##################### Description of haplotype combinations ########################

## exclude unclassified haplotypes of sam.AS, those not following Dan's rules

sam.ASE <- hap.d(sam.AS)

## table hap groups

sam.ASE[, ACDEF.pair:=paste0(ACDEF.1,"|",ACDEF.2)]

sam.ASE[,.N,by=ACDEF.pair]

##table of A or A- =1

sA <- sam.ASE[ACDEF.1 %in% grep("^1",sam.ASE$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE$ACDEF.2,value=T),.N, by=ACDEF.pair]

sum(sA$N)

## save and read

write.table(sam.ASE, file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/samples.ASE.txt", row.names=F)

sam.ASE <- fread(file="/mrc-bsu/scratch/ev250/ase-il2ra/objects/samples.ASE.txt")

########### strategy 2A
sam.ASE2 <- hap.d(sam.AS2)

sam.ASE2[, ACDEF.pair:=paste0(ACDEF.1,"|",ACDEF.2)]
##table of A or A- =1

sA2 <- sam.ASE2[ACDEF.1 %in% grep("^1",sam.ASE2$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE2$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE2$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE2$ACDEF.2,value=T),.N, by=ACDEF.pair]




######## strategy 2B

sam.ASE3 <- hap.d(sam.AS3)

sam.ASE3[, ACDEF.pair:=paste0(ACDEF.1,"|",ACDEF.2)]
##table of A or A- =1

sA3 <- sam.ASE3[ACDEF.1 %in% grep("^1",sam.ASE3$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE3$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE3$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE3$ACDEF.2,value=T),.N, by=ACDEF.pair]


###############################################################################
## plot ASE-Aprotective hap vs A hap combinations excluding samples with unclassified haplotypes (sam.ASE)

## extract haps starting with 1 or -

sam.A <-  sam.ASE[ACDEF.1 %in% grep("^1",sam.ASE$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE$ACDEF.2,value=T),]

## recode ASE%, from the hap.ASE col starting with 1 or -, if both choose 1

sam.A <- ASE.r(sam.A)

## PLot

## Add info on snps for phasing
wb <- loadWorkbook("/mrc-bsu/scratch/ev250/ase-il2ra/objects/CM_ASE_genotypes.xlsx")

##snps with tag only genotyping
tag.only <-data.table(readWorksheet(wb,sheet="All ichip SNPs",startRow=63,endRow=74,startCol=1, endCol=1, header=F))

sam.A[,Phasing:="ichip"][Donor %in% toupper(tag.only$Col1),Phasing:="Tag_SNPs"]

## get ymin and ymax for each hap pair group to add to plot

sam.A[,ACDEF.pair:=factor(ACDEF.pair, levels=unique(sam.A$ACDEF.pair))]
ymin <- sam.A[,min(ASE.hap2), by=ACDEF.pair]
ymax <- sam.A[,max(ASE.hap2), by=ACDEF.pair]
x <- sam.A[,mean(ASE.hap2),by=ACDEF.pair]

ggplot(sam.A,aes(x=ACDEF.pair, y=ASE.hap2, fill=Phasing, label=Donor)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))  + geom_hline(yintercept = 50,linetype="dashed") + geom_text(aes(label=ifelse(ASE.hap2<45, as.character(Donor),'')),hjust=1.2,vjust=1) + xlab("ACDEF(hap1|hap2)") + ylab("Hap2 ASE (%)") + theme(panel.grid.major.x = element_line(color = "grey80"))




## error bars incompatible with aes(fill)
## + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="red", width=0.2)

## check raw data for discordant haplotypes

sam.A[ACDEF.pair=="00110|10010",]

sam <- which(samples$V1 %in% sam.A[ACDEF.pair=="00110|10010",Donor])

col.hap.samp <- 6 + 2*(sam-1) ## get sample first hap in haps.samp

haps.samp[,sort(c(1:5,col.hap.samp, col.hap.samp+1)),with=F] ## get phasing data from shapeit output, seems OK

## compare with ped input

ped[ID %in% c(samples$V1[sam], "s10207617D"), c("ID",paste0(sort(rep(haps.D$SNP,2)), c(".1",".2"))), with=F]

####################### plot strategy 2A

sam.A2 <-  sam.ASE2[ACDEF.1 %in% grep("^1",sam.ASE2$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE2$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE2$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE2$ACDEF.2,value=T),]

## recode ASE%, from the hap.ASE col starting with 1 or -, if both choose 1

sam.A2 <- ASE.r(sam.A2)

## PLot

## Add info on snps for phasing
wb <- loadWorkbook("/mrc-bsu/scratch/ev250/ase-il2ra/objects/CM_ASE_genotypes.xlsx")

##snps with tag only genotyping
tag.only <-data.table(readWorksheet(wb,sheet="All ichip SNPs",startRow=63,endRow=74,startCol=1, endCol=1, header=F))

sam.A2[,Phasing:="ichip"][Donor %in% toupper(tag.only$Col1),Phasing:="Tag_SNPs"]

sam.A2[,ACDEF.pair:=factor(ACDEF.pair, levels=unique(sam.A2$ACDEF.pair))]


ggplot(sam.A2,aes(x=ACDEF.pair, y=ASE.hap2, fill=Phasing, label=Donor)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))  + geom_hline(yintercept = 50,linetype="dashed") + geom_text(aes(label=ifelse(ASE.hap2<45, as.character(Donor),'')),hjust=1.2,vjust=1) + xlab("ACDEF(hap1|hap2)") + ylab("Hap2 ASE (%)")  + theme(panel.grid.major.x = element_line(color = "grey80"))



############## plot strategy 2B #############

## extract haps starting with 1 or -

sam.A3 <-  sam.ASE3[ACDEF.1 %in% grep("^1",sam.ASE3$ACDEF.1,value=T) | ACDEF.1 %in% grep("^-", sam.ASE3$ACDEF.1,value=T) | ACDEF.2 %in% grep("^1",sam.ASE3$ACDEF.2,value=T) | ACDEF.2 %in% grep("^-", sam.ASE3$ACDEF.2,value=T),]

## recode ASE%, from the hap.ASE col starting with 1 or -, if both choose 1

sam.A3 <- ASE.r(sam.A3)

## PLot

## Add info on snps for phasing

wb <- loadWorkbook("/mrc-bsu/scratch/ev250/ase-il2ra/objects/CM_ASE_genotypes.xlsx")

##snps with tag only genotyping
tag.only <-data.table(readWorksheet(wb,sheet="All ichip SNPs",startRow=63,endRow=74,startCol=1, endCol=1, header=F))

sam.A3[,Phasing:="ichip"][Donor %in% toupper(tag.only$Col1),Phasing:="Tag_SNPs"]


sam.A3[,ACDEF.pair:=factor(ACDEF.pair, levels=unique(sam.A3$ACDEF.pair))]


ggplot(sam.A3,aes(x=ACDEF.pair, y=ASE.hap2, fill=Phasing, label=Donor)) + geom_dotplot(binaxis = "y", stackdir = "center", dotsize=0.5)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))  + geom_hline(yintercept = 50,linetype="dashed") + geom_text(aes(label=ifelse(ASE.hap2<45, as.character(Donor),'')),hjust=1.2,vjust=1) + xlab("ACDEF(hap1|hap2)") + ylab("Hap2 ASE (%)")  + theme(panel.grid.major.x = element_line(color = "grey80"))


## remove hap.pair and save sam.A3
sam.A3[,hap.pair:=NULL]

write.table(sam.A3, '/mrc-bsu/scratch/ev250/ase-il2ra/objects/haps.phasing.439snps.Ap495p.txt')
####################################################################################################
################################## Modelling effect of hapA as random effect ######################

## use sam.AS3 as input, most reliable phasing
## exclude samples with tag only genotyping

sam.AS3.e <- sam.AS3[!Donor %in% toupper(tag.only$Col1),]

## order samples so those starting with 1 or - are in the same col

sam.AS3.e <- hap.or(sam.AS3.e)

## recode samples: add new col A=0 if both haps are equal or - and 0. A=1 if 0 and 1 or - and 1

samAS3.r <- r.re(x=sam.AS3.e, cols=paste0("ACDEF", c(".1", ".2")))

## create grouping variable CDEF.g for running regression

samAS3.r[,CDEF.g:=paste0(CDEF.1,"|",  CDEF.2)]

## run regression, random intercept grouping by CDEF.g

re.lm <- lmer(ASE.hap2 ~ A + (1|CDEF.g), data = samAS3.r) 

write.table(samAS3.r, '/mrc-bsu/scratch/ev250/ase-il2ra/objects/haps.phasing.ran.eff.input.txt')

##plot
ggplot(samAS3.r,aes(x=CDEF.g, y=ASE.hap2, color=as.factor(A))) + geom_point()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))  + geom_hline(yintercept = 50,linetype="dashed") +  xlab("CDEF(hap1|hap2)") + ylab("Hap2 ASE (%)")  + theme(panel.grid.major.x = element_line(color = "grey80"))
