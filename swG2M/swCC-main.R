rm(list = ls())

library("ggplot2")
library("reshape2")
library("deSolve")
library("rootSolve")


setwd("C:/Users/Rosa/Downloads")

path <- "C:/Users/Rosa/Downloads"
spath <- "C:/Users/Rosa/Downloads"

source(paste(path, "fanalysis.R", sep = "/"))
source(paste(path, "swCC-models.r", sep = "/"))

source(paste(path, "swCC-general.R", sep = "/"))

mod = swcc.mod[[mdl]]
modPar = swcc.par[[pmdl]]
ic <- modPar$ic #inicial conditions of state variables
pp <- modPar$pp #initial parameters


ics <- do.call(cbind,lapply(molec, function(x,m,i){
  m <- m[,grep(x,colnames(m))]
  if(i[colnames(m)[1]]>=i[colnames(m)[2]]){
    m[,colnames(m)[1]] <- rev(st)
    m[,colnames(m)[2]] <- st
  }else{
    m[,colnames(m)[2]] <- rev(st)
    m[,colnames(m)[1]] <- st
  }
  return(as.matrix(m))
},ics,ic))


#bifurcation analysis
system.time(
  lsps <- ssa1p(initCon=ics,parSet=as.data.frame(t(pp)),mxInCon=st,val1par=vv1,nam1par=vn1, fmod=mod, tt=0, ss=0)
)

#save results
ba1p <- list(lsps, as.data.frame(t(pp)))
names(ba1p) <- c("lvalues","pps")
save(ba1p, file = paste(spath, sfile, sep="/"))

#######
#graphical representation
#######

#shape
gSPS <- melt(do.call(rbind,lsps), id.vars=c(vn1, "st"))
colnames(gSPS) <- c("par", "st", "variable", "value")
gSPS <- gSPS[which(gSPS$par>=0&gSPS$par<=30),]

head(gSPS)

#explore the stability
ggplot()+ theme_bw(base_size = 15)+
  geom_point(data=gSPS,aes(par,value, colour=factor(st)),size = 1.5)+#[which(gSPS$st%in%c(2,30)),]
  facet_wrap(~variable, nrow = 4)


#360x250
# plot the results we want
ggplot()+ theme_bw(base_size = 15)+
  theme(legend.title=element_blank()) +#+ scale_colour_manual(values = c("x"="#ff6600"))+
  xlab("[Chp] (AU)") +
  ylab("Steady State Concentration (AU)") +
  ggtitle("chk-G2M system")+
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%mVariables&gSPS$value<=0.5),], aes(par,value, colour=factor(variable)),size = 1.5) +
  geom_line(data=gSPS[which(gSPS$st%in%c(2,4)&gSPS$variable%in%mVariables&gSPS$value>=0.5),], aes(par,value, colour=factor(variable)),size = 1.5) +
  geom_line(data=gSPS[which(gSPS$st%in%c(0,3)&gSPS$variable%in%mVariables),], aes(par,value, colour=factor(variable)),size = 0.1) +#&gSPS$nt>=0.1&gSPS$nt<=6.07
  scale_y_continuous(breaks=seq(min(gSPS$value),pp["Z"],1), limit=c(0, pp["Z"])) +
  scale_x_continuous(breaks=seq(0,max(gSPS$par),5)) +
  scale_colour_manual(values=cols[which(names(cols)%in%mVariables)], 
                      breaks=mVariables,
                      labels=names(molec)[which(paste(molec, "0", sep = "")%in%mVariables)])
