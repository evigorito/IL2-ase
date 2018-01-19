library(data.table)
library(gtools)
#' match haplotype group
#'
#' This function allows you to match the haplotype type for a given input
#' @param haps DT with columns named as hap (haplotype group), SNP, protective and susceptible alleles
#' @param DT data table with phased information for a set of snps. rows individuals, 2 columns for each snp, corresponding to genotype in allele 1 and 2.
#' @keywords assign haplotype group
#' @export
#' @return input DT with new columns with hap info
#' m.hap()
m.hap <- function(haps,DT){
    hap.t <- unique(haps$hap)
     ## add cols for each allele
    DT[,sort(paste0(rep(hap.t,2), c(".1",".2"))):=lapply(sort(rep(hap.t,2)), function(i) rep("NA", nrow(DT)))] 
    ## identify missing values, coded 0
    m.v <- apply(DT[,  grep("rs",names(DT)), with=F], 2, function(i) which(i==0))
     if(length(m.v)>0){ # missing values true
             m.v <- m.v[unlist(lapply(m.v, function(i) length(i)>0))] # get snps and individuals
             m.snps <- unique(gsub("\\.1","",names(m.v))) # get unique missing snps ID
             m.g <- haps[SNP %in% m.snps,hap] # get hap group           
        }
        
    for(j in hap.t){ # need to do for allele 1 and allele 2
         snps <- haps[hap==j,SNP]
        r1 <- which(rowSums(DT[,snps,with=F] ==matrix(rep(haps[hap==j,protective],nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to protective allele1
        DT[r1,(paste0(j,".1")):=paste0(j,".p")]

        s1 <- which(rowSums(DT[,snps,with=F] ==matrix(rep(haps[hap==j,susceptible],nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to susceptible allele1
        DT[s1,(paste0(j,".1")):=paste0(j,".s")]

         r2 <- which(rowSums(DT[,paste0(snps, ".1"),with=F] ==matrix(rep(haps[hap==j,protective],nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to protective allele2
        DT[r2,(paste0(j,".2")):=paste0(j,".p")]

        s2 <- which(rowSums(DT[,paste0(snps, ".1"), with=F] ==matrix(rep(haps[hap==j,susceptible],nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to susceptible allele2
        DT[s2,(paste0(j,".2")):=paste0(j,".s")]
      
       
        if(j=="A"){ # define A-660.p/s (swaps)
           
            sp1 <- which(rowSums(DT[,snps,with=F] ==matrix(rep(c(haps[SNP==snps[1],susceptible], haps[SNP==snps[2],protective]),nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to protective allele1
            DT[sp1,A.1:="A-660.p"]
              
             ss1 <- which(rowSums(DT[,snps,with=F] ==matrix(rep(c(haps[SNP==snps[1],protective], haps[SNP==snps[2],susceptible]),nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to sus allele1
            DT[ss1,A.1:="A-660.s"]   

             sp2 <- which(rowSums(DT[,paste0(snps, ".1"), with=F] ==matrix(rep(c(haps[SNP==snps[1],susceptible], haps[SNP==snps[2],protective]),nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to protective allele1
            DT[sp2,A.2:="A-660.p"]
              
             ss2 <- which(rowSums(DT[,paste0(snps, ".1"),with=F] ==matrix(rep(c(haps[SNP==snps[1],protective], haps[SNP==snps[2],susceptible]),nrow(DT)), nrow=nrow(DT),byrow=T))==nrow(haps[hap==j,])) # matches to sus allele1
            DT[ss2,A.2:="A-660.s"]   
        }
                
        ##deal with missing values in D SNP rs56382813 coded as 0, extended to other groups
        if(exists("m.g")){
            if(j %in% m.g & !j %in% c("A","E","F")){ # exclude groups with one SNP or group A
                snps.nm <- snps[-which(snps %in% m.snps)]
                ## recode individuals in m.v
                ind <- m.v[[which(names(m.v) %in% haps[hap==j,SNP])]]
                
                mr1 <- which(rowSums(DT[ind,snps.nm,with=F] ==matrix(rep(haps[hap==j &SNP %in%  snps.nm ,protective],length(ind)), nrow=length(ind),byrow=T))==nrow(haps[hap==j & SNP  %in%  snps.nm,])) # matches to protective allele1
                DT[ind[mr1],(paste0(j,".1")):=paste0(j,".p")]

                ms1 <- which(rowSums(DT[ind,snps.nm,with=F] ==matrix(rep(haps[hap==j & SNP  %in%  snps.nm,susceptible],length(ind)), nrow=length(ind),byrow=T))==nrow(haps[hap==j& SNP  %in%  snps.nm ,])) # matches to susceptible allele1
                DT[ind[ms1],(paste0(j,".1")):=paste0(j,".s")]

                mr2 <- which(rowSums(DT[ind,paste0(snps.nm, ".1"),with=F] ==matrix(rep(haps[hap==j & SNP  %in%  snps.nm,protective],length(ind)), nrow=length(ind),byrow=T))==nrow(haps[hap==j & SNP  %in%  snps.nm,])) # matches to protective allele2
                DT[ind[mr2],(paste0(j,".2")):=paste0(j,".p")]

                ms2 <- which(rowSums(DT[ind,paste0(snps.nm, ".1"), with=F] ==matrix(rep(haps[hap==j & SNP  %in%  snps.nm,susceptible],length(ind)), nrow=length(ind),byrow=T))==nrow(haps[hap==j & SNP  %in%  snps.nm,])) # matches to susceptible allele2
                DT[ind[ms2],(paste0(j,".2")):=paste0(j,".s")]
            }
        }
    }
    
                    
    return(DT)
}


#' recode haplotype group
#'
#' This function allows you to recode the haplotype group:  1 for heterozygous; 0 for anything else
#' @param  DT with columns named ID and the columns created by m.hap (for the output of m.hap, cols<- unlist(lapply(unique(haps$hap), function(i) grep(paste0("^",i), x=names(DT), value=T)))), can be also full output from m.hap
#' @param groups hap groups to recode
#' @param DT with columns named as hap (haplotype group), SNP, protective and susceptible alleles
#' @keywords recode haplotype group
#' @export
#' @return data table with ID and recoded haps
#' r.hap()
r.hap <- function(DT,groups,haps){
    hap.cols <- unlist(lapply(unique(haps$hap), function(i) grep(paste0("^",i), x=names(DT), value=T)))
    tmp <- DT[,c("ID",hap.cols), with=F]
    tmp[,(groups):=0]
    if("A" %in% groups & "A-" %in% groups){
        tmp[A.1!=A.2 & A.1!="A-660.p" & A.2 !="A-660.p"& A.1!="A-660.s" & A.2 !="A-660.s", A:=1]
        tmp[A.1=="A-660.p" & A.2 =="A.s", `A-`:=1]
    }
    
    for(i in groups[which(!groups %in% c("A","A-"))]){
        cols <- grep(i,hap.cols, value=T)
        tmp[get(cols[1])==paste0(i,".p") & get(cols[2])==paste0(i,".s"), (i):=1] # p|s=1
        tmp[get(cols[1])==paste0(i,".s") & get(cols[2])==paste0(i,".p"), (i):=-1] #s|p=-1
    }
   
    
    return (tmp)
}

   
#' add haplotype info to sample file and recode using Dan's rules after phasing with shapeit
#'
#' This function allows you to add phasing info for each individual after phasing
#' @param  samples sample file output from phasing excluding first 2 rows
#' @param haps hap output from shapeit with SNP phasing info
#' @param haps.D data table with SNPs defining each hap and indicating protective or susceptible allele, column protec coded as 0 for ref allele and 1 for alt allele
#' @keywords phasing formatting
#' @export
#' @return data table with sample ID and phased SNPs recoded 1 for protective and 0 for susceptible (note that E hap susceptible is the minor allele)
#' phased.r()

phased.r <- function(samples, haps.samp,haps.D){
    s <- seq(6,ncol(haps.samp),2)
    h1 <- t(haps.samp[,paste0("V",s),with=F])
    h2 <- t(haps.samp[,paste0("V",s+1),with=F])
    colnames(h1) <- colnames(h2) <- haps.samp$V2
    ## order snps in h1 and h2 as in haps.D$SNP
    h1 <- data.table(h1[,colnames(h1)[order(match(colnames(h1),haps.D$SNP))]])
    h2 <- data.table(h2[,colnames(h2)[order(match(colnames(h2),haps.D$SNP))]])
    ## recode haps.D$protective as 0 or 1 based ref alt from haps.samp
    hap.D <- merge(haps.D, haps.samp[V2 %in% haps.D$SNP,.(V2,V4,V5)], by.x="SNP", by.y="V2", sort=FALSE)
    names(hap.D)[(ncol(hap.D)-1):ncol(hap.D)] <- c("ref","alt")
    hap.D[, protec:=ifelse(protective==alt,1,0)]
    
    ## define haps
    u <- unique(haps.D$hap)
    for(i in u){
        if(i=="A"){
            h1[,A:="0"][rs61839660==hap.D[SNP=="rs61839660",protec],A:="1"][rs61839660!=hap.D[SNP=="rs61839660",protec] & rs12722495==hap.D[SNP=="rs12722495",protec],A:="-"]
            h2[,A:="0"][rs61839660==hap.D[SNP=="rs61839660",protec],A:="1"][rs61839660!=hap.D[SNP=="rs61839660",protec] & rs12722495==hap.D[SNP=="rs12722495",protec],A:="-"]
        } else {
            p <- hap.D[hap==i,.(SNP,protec)]
            h1sub <- h1[,p$SNP,with=F]
            ## get whether h1sub has protective, susceptible or mixed haplotype
            hap1 <- apply(h1sub,1, function(j) sum(j==p$protec)/nrow(p)) ##0 is susc, 1 is protec,fraction is mixture 
            hap1[which(hap1>0 & hap1<1)] <- 9
            h1[,(i):=hap1]
            
            h2sub <- h2[,p$SNP,with=F]
            ## get whether h2sub has protective, susceptible or mixed haplotype
            hap2 <- apply(h2sub,1, function(j) sum(j==p$protec)/nrow(p)) ##0 is susc, 1 is protec,fraction is mixture 
            hap2[which(hap2>0 & hap2<1)] <- 9
            h2[,(i):=hap2]
        }
    } 
      
    u.hap <- paste0(u[-length(u)],collapse="")
    u.ASE <- paste0(u.hap,".",u[length(u)])
    h1[,(u.hap):=apply(h1[,(.SD), .SDcols = u[-length(u)]], 1, paste,collapse="")]
    h2[,(u.hap):=apply(h2[,(.SD), .SDcols = u[-length(u)]], 1, paste,collapse="")]

    h1[,(u.ASE):=apply(h1[,(.SD), .SDcols = u], 1, paste,collapse="")]
    h2[,(u.ASE):=apply(h2[,(.SD), .SDcols = u], 1, paste,collapse="")]
    
    ## add haps to samples
    samples[,(paste0(u.hap,".1")):=h1[[u.hap]]][,(paste0(u.hap,".2")):=h2[[u.hap]]]
    samples[,(paste0(u.ASE,".1")):=h1[[u.ASE]]][,(paste0(u.ASE,".2")):=h2[[u.ASE]]]
    
    
    return(samples)
}

#' format haps file (output from shapeit) as ped like
#'
#' This function allows you to reformat haps output file from shapeit as ped like
#' @param  samples sample file output from phasing excluding first 2 rows, or any file with first col named V1 with sample order
#' @param hap hap output from shapeit with SNP phasing info
#' @keywords haps formatting
#' @export
#' @return data table with sample ID and snps info as ped
#' haps.ped()

haps.ped <- function(samples,hap){
    s <- seq(6,ncol(hap),2)
    h1 <- data.table(t(hap[,paste0("V",s),with=F]))
    h1 <- h1[,lapply(.SD, as.character)]
    h2 <- data.table(t(hap[,paste0("V",s+1),with=F]))
    h2 <- h2[,lapply(.SD, as.character)]
    ##recode
    DT <- hap[,c("V2","V4","V5"),with=F]
    ## for(i in 1:nrow(DT)){
    ##     h1[get(paste0("V",i))==0,paste0("V",i):=DT[i,V4]][get(paste0("V",i))==1,paste0("V",i):=DT[i,V5]]
    ##     h2[get(paste0("V",i))==0,paste0("V",i):=DT[i,V4]][get(paste0("V",i))==1,paste0("V",i):=DT[i,V5]]
    ## }
    ## name h1 and h2
    names(h1) <- paste0(DT$V2,".1")
    names(h2) <- paste0(DT$V2,".2")
    ## combine h1 and h2 and add sample ID
    tmp <- cbind(samples$V1,h1,h2)
    ##order names tmp
    setcolorder(tmp,names(tmp)[c(1,as.vector(sapply(DT$V2, function(i) grep(i,names(tmp)))))])
    names(tmp)[1] <- "ID"
    return(tmp)
}

#' recode ped file with genotypes 0,1, NA
#'
#' This function allows you to recode ped genotypes as numeric
#' @param  ped1 ped file to recode
#' @param hap data table with V2 col= rsid, V4= ref allele and V5=alt allele, can be hap file
#' @keywords ped recoding
#' @export
#' @return data table with sample ID and snps coded 0,1
#' r.ped()

r.ped <- function(ped1,hap){
    DT <- copy(ped1) ## so ped1 is not modified
    for(i in 1:nrow(hap)){
        g <- grep(hap[i,V2],names(DT))
        for(j in 1:2){ ## two cols per snp
            DT[get(names(DT)[g[j]])==0,names(DT)[g[j]]:="NA"][get(names(DT)[g[j]])==hap[i,V4],names(DT)[g[j]]:="0"][get(names(DT)[g[j]])==hap[i,V5],names(DT)[g[j]]:="1"]
        }
    }
    return(DT)
}


#' convert ped/map into GEN/sample format
#'
#' This function allows you to convert ped/map into GEN/sample
#' @param ped ped file
#' @param map map file
#' @keywords ped GEN conversion
#' @export
#' @return list with datatables in GEN and sample formats
#' ped2gen()

ped2gen <- function(ped,map){
    gen <- map[,.(Chr,refsnp_id,b37,a0,a1)]
    g <- unique(gsub("\\.[1-2]$","",names(ped)[7:ncol(ped)]))
    for(i in 1:nrow(gen)){
        p <- permutations(n=2,r=2,v=unlist(gen[i, .(a0,a1)]),repeats.allowed=T)
        ################ work in progress #########
    }
}


      
#' compare ped with hap file after recoding
#'
#' This function allows you to recode ped genotypes as numeric
#' @param  ped1 ped file recoded to 0,1,NA output from r.ped
#' @param h.ped hap file in ped format output from haps.ped
#' @param snps vector with snps to compare
#' @keywords ped hap comparison
#' @export
#' @return names vector with snp list, each entry is the proportion of consistency: 1 is 100% consistency. NA when one of the files has NAs (input for phasing)
#' h.p.comp()

h.p.comp <- function(ped1,h.ped,snps){
    v <- c()
    for (i in snps.d){
        ped1 <- ped1[, lapply(.SD,as.numeric)]
        h.p <- h.ped[,2:ncol(h.ped), with=F][,lapply(.SD,as.numeric)]
        ped.s <- rowSums(ped1[, grep(i, names(ped1),value=T),with=F])
        h.s <- rowSums(h.p[, grep(i, names(h.p),value=T),with=F])
        v <- c(v,sum(ped.s==h.s)/length(ped.s))
    }
    names(v) <- snps.d
    return(v)
}

#' order haps so those starting with 1 or - are on the same column
#'
#' This function allows you to order haps so those starting with 1 or - are on the same column
#' @param  sam.ASE output from phased.r with ASE col and haps starting with 1 or -
#' @keywords hap order
#' @export
#' @return input DT with ordered haps
#' hap.or()

hap.or <- function(sam.ASE){
    pat <- c("1","-")
    g <- grep("ASE",names(sam.ASE),value=T) # cols to swap with 
    w <- lapply(g, function(i) sapply(pat, function(j) grep(paste0("^",j), sam.ASE[[i]])))
    ## get rows to swap
    sw <- c(w[[1]]$`1`[which(!w[[1]]$`1` %in% w[[2]]$`1`)], w[[1]]$`-`[which(!w[[1]]$`-` %in% unlist(w[[2]]))]) 
    to.sw <- sam.ASE[sw,] ## keep values to ease swapping
    ## identify cols to swap
    pat2 <- c(".1",".2")
    scol <- sapply(pat2, function(i) grep(paste0(i,"$"),names(sam.ASE), value=T))
    sam.ASE[sw, scol[,1]:=lapply(1:nrow(scol), function(i) to.sw[[scol[i,2]]])]
    sam.ASE[sw, scol[,2]:=lapply(1:nrow(scol), function(i) to.sw[[scol[i,1]]])]
    setkeyv(sam.ASE,scol[1,2:1])
    return(sam.ASE)

}



#' recode ASE%, assign ASE (rs.ASE,ref,alt col) to hap2
#'
#' This function allows you to assign ASE % to hap2 (haps previously swapped using hap.or)
#' @param  sam.A output from phased.r and hap.or 
#' @keywords ASE recode
#' @export
#' @return input DT with an extra column with recoded ASE
#' ASE.r()

ASE.r <- function(sam.A){
    sam.A[,ASE.hap2:=rs12244380.G.A]
    ## recode for hap2.ASE finishing in 1
    g <- grep("ASE.2",names(sam.A),value=T)
    rec <- grep("1$", sam.A[[g]])
    sam.A[rec,ASE.hap2:=100-ASE.hap2]
    return(sam.A)
}

    
#' Description of haplotype combinations after phasing
#'
#' This function allows to do format haplotype combination to ease plotting and presentation in tables
#' @param sam.AS data table with phasing and ASE data after manipulation of phasing output 
#' @keywords ASE hap descritpion
#' @export
#' @return DT with ordered hap combinations. Excludes unclassified haps (dont follow Dans rules)
#' hap.d()

hap.d <- function(sam.AS){
    
    mis <- unique(as.vector(apply(sam.AS[,(ncol(sam.AS)-3):ncol(sam.AS),with=F], 2, grep, pattern="9"))) ## remove cols with 9

    sam.ASE <- sam.AS[!mis,]
    
    ## order haps so those starting with 1 or - are on the same column
    
    sam.ASE <- hap.or(sam.ASE)

    ## reformat

    sam.ASE <- hap.sub(sam.ASE)
       
    return(sam.ASE)
}

#' Subfunction for description of haplotype combinations after phasing
#'
#' This function allows to do format haplotype combination to ease plotting and presentation in tables
#' @param sam.ASE data table with phasing and ASE data after manipulation of phasing output
#' @keywords ASE hap descritpion
#' @export
#' @return DT with ordered hap combinations. Excludes unclassified haps (dont follow Dans rules)
#' hap.sub()

hap.sub <- function(sam.ASE){
    ## get the unique pairs of haplotypes from the set of haplotypes in sam.ASE

    uhaps <- unique(c(sam.ASE[[names(sam.ASE)[4]]], sam.ASE[[names(sam.ASE)[5]]])) ##unique haps
    u.pairs <- u.hap.pairs2(uhaps) ## unique hap pairs
    
    ##u.pairs[,Var1:=factor(Var1, levels= unique(sam.ASE$ACDEF.2)) ]
    
    ##u.pairs[,Var2:=factor(Var2,levels=  unique(sam.ASE$ACDEF.2))]
    
    ##setkey(u.pairs,Var2,Var1)

    ## for each individual in sam.ASE get the u.pairs row number that matches its hap.pair
    sam.ASE[,hap.pair:=0]
    for(i in 1:nrow(sam.ASE)){
        w1 <- which( u.pairs$Var1  %in% sam.ASE[i,ACDEF.1])
        w2 <- which(u.pairs$Var2  %in% sam.ASE[i,ACDEF.2])
        w <- which(w1 %in% w2)
        if(length(w)!=0){
            sam.ASE[i,hap.pair:=w1[w]]
        } else { #can be swapped
            w1 <- which(u.pairs$Var2  %in% sam.ASE[i,ACDEF.1])
            w2 <- which(u.pairs$Var1  %in% sam.ASE[i,ACDEF.2])
            w <- which(w1 %in% w2)
            sam.ASE[i,hap.pair:=w1[w]]
            ## swap to ease identifying haps pairs
            s <- sam.ASE[i,ACDEF.1]
            s2 <- sam.ASE[i,ACDEF.ASE.1]
            sam.ASE[i, ACDEF.1:=ACDEF.2][i,ACDEF.2:=s]
            sam.ASE[i,ACDEF.ASE.1:=ACDEF.ASE.2][i,ACDEF.ASE.2:=s2]
        }
    }

    setkey(sam.ASE,hap.pair)
    
    return(sam.ASE)
}

#' Subfunction for description of haplotype combinations after phasing, different input
#'
#' This function allows to do format haplotype combination to ease plotting and presentation in tables
#' @param sam.ASE data table with phasing and ASE data after manipulation of phasing output
#' @param cols name of cols to look for hap pairs
#' @keywords ASE hap description
#' @export
#' @return DT with ordered hap combinations. Excludes unclassified haps (dont follow Dans rules)
#' hap.sub2()

hap.sub2 <- function(sam.ASE,cols=NULL){
    ## get the unique pairs of haplotypes from the set of haplotypes in sam.ASE
    if(is.null(cols)){
        cols <- names(sam.ASE)[4:5]
    } 
    uhaps <- unique(c(sam.ASE[[which(names(sam.ASE)==cols[1])]], sam.ASE[[which(names(sam.ASE)==cols[2])]]))    
    u.pairs <- u.hap.pairs2(uhaps) ## unique hap pairs
    
    ##u.pairs[,Var1:=factor(Var1, levels= unique(sam.ASE$ACDEF.2)) ]
    
    ##u.pairs[,Var2:=factor(Var2,levels=  unique(sam.ASE$ACDEF.2))]
    
    ##setkey(u.pairs,Var2,Var1)

    ## for each individual in sam.ASE get the u.pairs row number that matches its hap.pair
    sam.ASE[,hap.pair:=0]
    for(i in 1:nrow(sam.ASE)){
        w1 <- which( u.pairs$Var1  %in% sam.ASE[i,cols[1], with=F])
        w2 <- which(u.pairs$Var2  %in% sam.ASE[i,cols[2], with=F])
        w <- which(w1 %in% w2)
        if(length(w)!=0){
            sam.ASE[i,hap.pair:=w1[w]]
        } else { #may be swapped
            w1 <- which(u.pairs$Var2  %in% sam.ASE[i,cols[1],with=F])
            w2 <- which(u.pairs$Var1  %in% sam.ASE[i,cols[2],with=F])
            w <- which(w1 %in% w2)
            sam.ASE[i,hap.pair:=w1[w]]
            ## swap to ease identifying haps pairs
            ## get ASE cols
            ase.cols <- grep("ASE",names(sam.ASE),value=T)
            s <- sam.ASE[i,cols[1],with=F]
            s2 <- sam.ASE[i,ase.cols[1],with=F]

            ## working here, recode swap by reference
            sam.ASE[i, cols[1]:=get(cols[2])][i,cols[2]:=s]
            sam.ASE[i,ase.cols[1]:=get(ase.cols[2])][i,ase.cols[2]:=s2]
        }
    }

    setkey(sam.ASE,hap.pair)
    
    return(sam.ASE)
}


#' Recode haps for modelling random effects: isolate A effect
#'
#' This function allows to recode haps as A plus CDEF
#' @param x data table with phasing and ASE data after manipulation of phasing output (sam.AS)
#' @param cols name of columns to split by A hap (i.e. "ACDEF.1", "ACDEF.2")
#' @keywords random effects hap recoding
#' @export
#' @return DT with recoded haps
#' r.re()

r.re <- function(x, cols){
    ## code haplotype A as 1 (het) or 0 (hom), "-" is equivalent to "0"
    x[,A:=0][(substr(ACDEF.1,1,1)==0 |substr(ACDEF.1,1,1)=="-") & substr(ACDEF.2,1,1)==1, A:=1]
    ## make new columns excluding A
    new.cols <- substr(cols,2,nchar(cols[1]))
    x[,(new.cols):= list(substr(get(cols[1]),2,nchar(cols[1])), substr(get(cols[2]),2,nchar(cols[1])))]
    ase.cols <- grep("ASE", names(x), value=T)
    new.ase <- substr(ase.cols,2,nchar(ase.cols[1]))
    x[,(new.ase):= list(substr(get(ase.cols[1]),2,nchar(ase.cols[1])), substr(get(ase.cols[2]),2,nchar(ase.cols[1])))]
    ## when A=0, swap new.cols to order the haplotypes
    tmp <- x[A==0, c("Donor", "Condition", "rs12244380.G.A", new.cols,new.ase), with=F]
    tmp2 <- hap.sub2(tmp)
    tmp2[,hap.pair:=NULL]
    ## add A col to tmp2
    tmp2[,A:=0]
    setcolorder(tmp2, c("Donor", "Condition", "rs12244380.G.A", "A", new.cols,new.ase))
    ## select c("Donor", "Condition", "rs12244380.G.A", new.cols,new.ase) from x and rbind x[A==1,] with tmp2 (excluding hap.pair)
    y <- rbind(tmp2,  x[A==1, c("Donor", "Condition", "rs12244380.G.A", "A", new.cols,new.ase), with=F])
    ## recode ASE in y relative to hap2
    z <- ASE.r(y)
    return(z)
}

    
