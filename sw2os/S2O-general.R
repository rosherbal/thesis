#main paths
path <- "C:/Users/Rosa/Downloads"
spath <- "C:/Users/Rosa/Downloads"

#saving file
sfile <- "s2o-bd1.Rdata"

#model
mdl <- "f5"
pmdl <- "bd1"

#time course analysis
tt <- c(0.0, 300,0.001)#t0, tfinal, step

#bifurcation analysis
mol <- "PP" #selected molecule to track
vn1 <- "nt" #names of parameter 1 to change
vv1 <- seq(0,60,0.1) #define parameter points to check #c(seq(0,10,0.02),seq(10.5,100,0.2))#c(seq(0,0.4,0.02),seq(0.5,2,0.1),seq(3,52,1),seq(52.05,53,0.05),seq(54,70,1),seq(70.5,100,0.5))#seq(0,10,0.05) # vv1 <- unique(ba1p[["lvalues"]][[1]]$nt) #seq(min(subset(ssp, cor == 2)$par),max(subset(ssp, cor == 1)$par),0.001)
jac <- T#calculate the jacobian matrix?
st <- seq(0.0,4,0.05) #set of initial condidions of state variables to check
sys <- F #proceed with a systematical analysis

#graphical conditions
labs <- c("Phosphate Donor (nt) (AU)", "Concentration (AU)", "")#c("Time (AU)", "Concentration (AU)", "")#xlab, ylab, title#
cols <- c("2" = "#619cff","0" ="#f8766d","3" = "#f564e3", "4" = "#00bfc4")#color scheme for types of point
cols2 <- c("OO" = "#ae0033", "PP" = "#006d8f", "OP" = "#9d9d9d", "PO" = "#5e5e5e")#color scheme for molecule#"OO" = "#ffa500",

#matrix of initial conditions from st
if (mdl == "f5"){nspecies <- c("OO","PO","PP")}else{nspecies <- c("OO","PP")}

ics <- matrix(data = NA, nrow = length(st), ncol = length(nspecies), dimnames = list(NULL,nspecies))
ics[,"OO"] <- st
ics[,"PP"] <- rev(st)
if (mdl == "f5"){
  ics[,"PO"] <- max(st) - (st + rev(st))  
}

#remove extream values min(st) and max(st)
ics[which(ics[,2]<0),2] <- 0
ics<-ics[-1,]
ics<-ics[-dim(ics)[1],]

#fix inicial conditions
ic <- ics[38,] 

#produce all combinations of parameters for systematical analysis

if(sys){
  pps <- expand.grid(0:1,0:1,0:1,0:1,0:1,0:1,0:1,0:1,c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05),c(0,0.05))
  colnames(pps) <-c("p0", "p1", "p2", "p3", "d0", "d1", "d2", "d3", "bp0", "bp1", "bp2", "bp3", "bd0", "bd1", "bd2", "bd3")
  mainp <-c("p0", "p1", "p2", "p3", "d0", "d1", "d2", "d3")
  
  if(pmdl=="ti"){

    # pps <- do.call(rbind,apply(pps,1,function(r){sr <- length(r)-length(r[which(r!=0)]);if(sr==2){return(r)}}))
    # pps <- t(sapply(mainp, function(p,dfpp){return(dfpp[which(dfpp[,p]==0 & dfpp[,paste("b",p, sep = "")]==0),])},pps))
    
     
  }# # single reaction path TI
  else if(pmdl=="ai"){
    
    pps <- do.call(rbind,apply(pps,1,function(r){sr <- length(r)-length(r[which(r!=0)]);if(sr==4){return(r)}}))
    pps <- (unique(do.call(rbind,lapply(mainp[grep("p", mainp)], function(p,dfpp,mainp){
      temp <- t(sapply(mainp[grep("d", mainp)],function(p2,p,dfpp){
        if(p2 != p){
          return(dfpp[which(dfpp[,p]==0 & dfpp[,paste("b",p, sep = "")]==0 & dfpp[,p2]==0 & dfpp[,paste("b",p2, sep = "")]==0),])}
        },p,pps))
      },dfpp,mainp))))
    
  }# # two opposite reaction path TI
  
  else if(pmdl=="bp"){
    pps <- pps[which(pps[,"p2"]==0 & pps[,"d0"]==0 &
                       pps[,"bp2"]==0 & pps[,"bd0"]==0),]#BP
    pps <- do.call(rbind,apply(pps,1,function(r){sr <- length(r)-length(r[which(r!=0)]);if(sr==5){return(r)}}))
  
  }# # one reaction from BP 
  else if (pmdl=="bd1"){
    pps <- pps[which(pps[,"p2"]==0 & pps[,"d0"]==0 &
                       pps[,"bp2"]==0 & pps[,"bd0"]==0 & pps[,"bd1"]==0),]#BD1
    pps <- do.call(rbind,apply(pps,1,function(r){sr <- length(r)-length(r[which(r!=0)]);if(sr==6){return(r)}}))
    
  }# # one reaction from BD1
  
}






