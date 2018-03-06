rm(list = ls())

library("ggplot2")
library("reshape2")
library("deSolve")
library("rootSolve")

setwd("C:/Linux/data/sw2os")

path <- "C:/Linux/data/sw2os"

source(paste(path, "fanalysis.R", sep = "/"))
source(paste(path, "S2O-general.R", sep = "/"))
source(paste(path, "S2O-models.R", sep = "/"))

#####
#variables
#####

mod = s2o.mod[[mdl]]
modPar = s2o.par[[pmdl]]
pp <- modPar$pp #initial parameters

#extend parameters of systematical analysis
if (sys){
  extpp <-cbind(rep(pp["nt"],dim(pps)[1]),rep(pp["nd"],dim(pps)[1]))
  extpp <- cbind(extpp,rep(max(st),dim(pps)[1]))
  colnames(extpp) <- c("nt","nd", "tot")
  pps <- as.data.frame(cbind(pps, extpp), stringsAsFactors = F)
  rm(extpp)
}else{
  pps <- t(data.frame(pp, stringsAsFactors = F))
}

#####
#### time course diagrams
##############

# pp["nt"]<- 15 #specific value for nucleotide donor
tc <- tcf(ic=ic,tt=tt,mod=mod,pp=pp)#funciton in fanalysis.R
tc$OP <- pp["tot"] - apply(tc[,-1],1, sum) #add the missing species
tc$OP[tc$OP<0.0] <- 0.0
gtc <- melt(as.data.frame(tc), id.vars="time")

#plot it!
ggplot(data=gtc,aes(time,value,colour=factor(variable)))+
  geom_line(size=1) +
  #geom_point(data=ver, aes(as.numeric(t),as.numeric(v.PO),colour=factor(n)), size=3)+
  scale_y_continuous(breaks=seq(0,pp["tot"],1), limit=c(0, pp["tot"])) +
  #scale_x_continuous(breaks=seq(min(tt),max(tt),50)) +
  scale_colour_manual(values = cols2, name="State")+
  xlab("Time (AU)") + ylab(labs[2]) + ggtitle(paste("nt = ",pp["nt"], sep = "")) + theme_bw(base_size = 20)


#### bifurcation analysis
##############

system.time(
  lsps <- ssa1p(initCon=ics,parSet=pps,mxInCon=st,val1par=vv1,nam1par=vn1, fmod=mod, ss=0)#lsps
)

#calculate bistable region for each paramSet
bsregion <- lapply(lsps,BiOsreg, mol=mol,nam1par=vn1, cutoffval=0.5)

#know if they are switches or oscillators
type <- do.call(c,lapply(lsps,function(df){cls <- unique(df$st); print(cls);
if(all(2%in%cls & 0%in%cls)){if(all(4%in%cls & 3%in%cls)){cond = "so"} else if(4%in%cls){cond = "ds"}else if(3%in%cls){cond = "hs"}else{cond = "s"}}
else if(3%in%cls){if(4%in%cls){cond = "do"}else{cond = "o"}}
else{cond=NA}; return(cond)}))

#save raw data
ba1p <- list(lsps,bsregion,type,pps)
names(ba1p) <- c("lvalues","lbsregion","type","pps")
save(ba1p, file = paste(spath, sfile, sep="/"))

#identify pps with removing parameters
outvals <- c("p2","d0","bd0","bp2")#c("p3","bp3")#
dict <- as.data.frame(cbind(1:length(ba1p[["lvalues"]]), apply(apply(ba1p[["pps"]],1,function(r,n){n[which(r==0)]},names(ba1p[["pps"]])),2,paste, collapse=" ")),stringsAsFactors = F)
dict$V3 <- sapply(dict[,2],function(j,outvals){j <-unlist(strsplit(j," ",fixed = T));k<-paste(j[-which(j%in%outvals)],collapse = " ");return(k)},outvals)
#calculate amplitude and period for oscillatory models (so,o,do types)

#select the matrixes
oscmod <- ba1p[["lvalues"]][which(ba1p[["type"]]%in%c("so","o","do"))]#grep("bd1", dict[,2])
names(oscmod)<-which(ba1p[["type"]]%in%c("so","o","do"))

if(length(oscmod) > 0){
  #calculate A and T for each set 
    AT3 <- lapply(i,function(i,pps,infr,tms,inc,thrs=5,period=F, mod=mod){#AT names(oscmod)#13#4.5#AT[i]
    i<- as.numeric(i)
    inf_up <- max(ba1p[["lvalues"]][[i]][which(ba1p[["lvalues"]][[i]]$st==3),"nt"]) + 0.2
    inf_dw <- min(ba1p[["lvalues"]][[i]][which(ba1p[["lvalues"]][[i]]$st==3),"nt"])- 0.2
    print("#############")
    print(i)
    print(c(inf_dw,inf_up))
    print("#############")

    if((inf_up-inf_dw)<=2){infr<- seq(inf_dw,inf_up,0.1)
    }else{infr <- seq(5,53,0.5)}#c(seq(inf_dw,(inf_dw+2), 0.1),seq((inf_dw+2.1), inf_up,0.5))}#seq(5.6,15,0.5)}#{
    # infr <- seq(inf_dw,inf_up,0.5)#0.5 
    if(min(infr)<0){infr<-infr[which(infr>0)]; infr<-c(0.01,infr)}
    m <- ATfun(inframe=infr,parset=pps[i,],timeseq=tms,incon=inc,thrs,period=F)
    m[,"gt"] <- i
    
   # save(m,file=paste("C:/Users/rosherbal/Dropbox/Kings/results/s2o/","ba1p_bp_up2one_AT_", i, ".Rdata", sep = "")) 
    
    return(m)},ba1p[["pps"]],ba1p[["bsregion"]],tt,ics[41,])

    names(AT)<- names(oscmod)
    
    save(AT,file=paste(spath,"ba1p_bp_up2one_AT.Rdata", sep = "/")) 

}

dfAT <- do.call(rbind,AT)

dfAT <- dfAT[which(dfAT$p<=15.1),]

for(i in unique(dfAT$gt)){dfAT$gt[which(dfAT$gt==i)] <- dict[i,3]}


#remove miscalculations for bistable models
lsps <- lapply(1:length(ba1p[["lvalues"]]), rmmiscal, lsps=ba1p[["lvalues"]], bsre=ba1p[["lbsregion"]], nam1par=vn1)
lsps <- lapply(1:length(lsps),function(n,l){cbind(l[[n]],gt=rep(n,dim(lsps[[n]])[1]))},lsps)


##prueba para 1 df
#ST1df <- lapply(3, rmmiscal, lsps=lsps, bsre=bsregion, nam1par=vn1)
#gST1 <- melt(ST1df, id.vars=c("nt", "st"))



#di to sp topologies
# dfsps$gt <-factor(dfsps$gt, levels=c("p0 d0 bp0 bd0", "p1 d1 bp1 bd1", "p2 d2 bp2 bd2", "p3 d3 bp3 bd3",
#                                      "p0 d1 bp0 bd1", "p1 d0 bp1 bd0", "p2 d3 bp2 bd3", "p3 d2 bp3 bd2",
#                                      "p0 d3 bp0 bd3", "p1 d2 bp1 bd2", "p2 d1 bp2 bd1", "p3 d0 bp3 bd0",
#                                      "p0 d2 bp0 bd2", "p1 d3 bp1 bd3", "p2 d0 bp2 bd0", "p3 d1 bp3 bd1"))

#preparing for ploting...

dfsps <- do.call(rbind,lsps)

for(i in unique(dfsps$gt)){
  dfsps$gt[which(dfsps$gt==i)] <- dict[i,3]
}

gSPS <- melt(dfsps, id.vars=c("nt", "st", "gt"))
gSPS<-gSPS[which(gSPS$nt<=15.1),]#15.1

#extra settings

#di_up2two
gSPS <- gSPS[which(gSPS$gt%in%c(dict[which(dict[,1]%in%c(16,15,14,13,8,7,6,5)),2])),]#di_up2two
#pps8_up2one
gSPS$gt <-factor(gSPS$gt, levels=c(dict[16,2],dict[15,2],dict[14,2],dict[13,2],
                                   dict[8,2],dict[7,2], dict[6,2],dict[5,2]))#pps8_up2one

# gSPS$gt <-factor(gSPS$gt, levels=c(dict[12,3],dict[11,3],dict[10,3],
#                                     dict[6,3],dict[5,3],dict[4,3],
#                                     dict[9,3],dict[8,3],dict[7,3],
#                                     dict[3,3], dict[2,3],dict[1,3]))#pps8_up2one
gSPS$gt <-factor(gSPS$gt, levels=c("p0","p1","p3",
                                   "bp0","bp1","bp3",
                                   "d1","d2","d3",
                                   "bd2","bd3"))#pps8_bd1_up2two

gSPS$gt <-factor(gSPS$gt, levels=c("p0","p1","p3",
                                   "bp0","bp1","bp3",
                                   "d1","d2","d3",
                                   "bd1","bd2","bd3"))#bp_up2one


# gSPS$gt <-factor(gSPS$gt, levels=c(dict[16,2],dict[15,2],dict[14,2],dict[13,2],
#                                    dict[8,2], dict[7,2],dict[5,2]))
# dfsps$gt <-factor(dfsps$gt, levels=c("p0 p2 d0 bp2 bd0 bd1","p1 p2 d0 bp2 bd0 bd1","p2 p3 d0 bp2 bd0 bd1",
#                                      "p2 d0 d1 bp2 bd0 bd1","p2 d0 d2 bp2 bd0 bd1","p2 d0 d3 bp2 bd0 bd1",
#                                      "p2 d0 bp0 bp2 bd0 bd1","p2 d0 bp1 bp2 bd0 bd1","p2 d0 bp2 bp3 bd0 bd1",
#                                      "p2 d0 bp2 bd0 bd1 bd2","p2 d0 bp2 bd0 bd1 bd3"))

##plot it

ggplot()+ #[gST1$st%in%c(2,4,3),]
  geom_line(data=gSPS[which(gSPS$st%in%c(3)&gSPS$variable%in%c("OO","PP")),], aes(nt,value, colour=factor(variable)), size = 1, linetype = "dashed" ) + #longdash, 2"dashed
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%c("OO","PP")),], aes(nt,value, colour=factor(variable)),size = 1.5) +#&gSPS$value>=0.5
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%c("OO","PP")&gSPS$value>=0.5),], aes(nt,value, colour=factor(variable)),size = 1.5) +
  geom_line(data=dfAT[which(dfAT$n%in%c("OO","PP")),], aes(as.numeric(p), as.numeric(An), colour=factor(n)), size = 0.5) +#linetype = "dotted"
  geom_line(data=dfAT[which(dfAT$n%in%c("OO","PP")),], aes(as.numeric(p),as.numeric(Ax), colour=factor(n)), size = 0.5) + #linetype = "dotted"
  facet_wrap(~gt, nrow = 4) +
  scale_y_continuous(breaks=seq(0,as.numeric(pp["tot"])), limit=c(0, as.numeric(pp["tot"]))) +
  scale_x_continuous(breaks=seq(0,max(gSPS$nt),5)) +
  scale_colour_manual(values = cols2, name="Form")+
  xlab(labs[1]) + ylab(labs[2]) + theme_bw(base_size = 15) #+ # + theme( panel.panel.spacing.y = unit(1, "lines")) #+ guides(colour=F)gST1 <- melt(ST1df, id.vars=c(vn1, "st"))#"par",
  ggtitle("BP system")