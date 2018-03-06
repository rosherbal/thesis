
#creat time course diagrams
# ic:initial condition set
# tt: time segment
# mod: function
# pp: parameter set
tcf <- function(ic,tt,mod,pp){
  tc <-  as.data.frame(ode(ic,seq(tt[1], tt[2], by = tt[3]), mod, pp,  maxsteps = 10000))
  #tc <-  as.data.frame(lsode(y=ic,times=seq(tt[1], tt[2], by = tt[3]), func=mod, params=pp,  rtol = 1e-8, maxsteps = 5000))
  return(tc) 
}

#steady state iteration
# x: parameter set
# nx: name of the parameter to change
# ic:initial condition set
# mod: function
# tt: time segment
sse <- function(x, nx, ics, mod, tt){
  temp <- as.data.frame(t(apply(ics,1,function(sv,x,mod,tt){#sv>initial conditions, x>parameters, mod>model
    st <- ode(y = sv, time = tt, func = mod, parms = x)#calculate the ODE
    st <- replace(st[nrow(st),-1],which(st[nrow(st),-1]<0),0)#Taking the last value
    return(c(st, x[which(names(x)%in%nx)], st = NA))#return the values
  },x,mod,tt)))#return the values to main
  return (temp)
}


#calculate the value at equilibrium through an interative solver. It clasify the equilibrium points through the calculation of the jacobian.
# x: parameter set
# nx: names of the parameter set
# ic:initial condition set
# mod: function
# tt: point of time (0 by default)
ssi <- function(x, nx, ics, mod, tt){
  s <- NULL
  temp <- as.data.frame(t(apply(ics,1,function(sv, x, mod, tt){ 
    #calculate the equilibrium point
    st <- withCallingHandlers(
      stode(y=sv, time = tt, func=mod, parms=x, maxiter = 1000, rtol = 1e-8, positive = T)$y,
      warning = function(w){
       if(grepl("steady-state not reached", w$message)| grepl("error during factorisation", w$message)){NULL#
          } else {message(w$message)}
      })
    #calculate the jacobian and clasify the equilibrium point accordingly
    if(is.null(st)){return(NA)}else{
     if(jac){
        jst <- eigen(jacobian.full(y = st, func = mod, parms = x))$values
        if(is.complex(jst)){
          jst <- sign(Re(jst))
          if (all(jst < 0)){s = 4}#stable
          else if(all(jst > 0)){s = 5}#unstable
          else if (any(jst < 0) && any(jst > 0)){s = 3}#HB
          else if (table(jst==0)["TRUE"]>=1){s = 6} 
        }else{#if it is not complex number (switches)
          jst <-  sign(jst)#sign of the eigenvalues
          if (all(jst < 0)){s = 2} #stable
          else if(all(jst > 0)){s = 1}#unstable
          else if (any(jst < 0) && any(jst > 0)){s = 0}#SN
          else if (table(jst==0)["TRUE"]>=1){s = 7} 
        }
      }else if(!jac){s = NA}
      #if (length(c(st, par = p, st = s))<5){}
      return(c(st, x[which(names(x)%in%nx)], st = s))
    }
  }, x, mod, tt)))
  return (temp)
}

#call for ssi or sse to calculate the equilibrium point (and for the ssi the jacobian). Create a df for each parSet and returns it.
# initCon: initial condition set
# parSet: parameter set
# mxInCon: range of values of initial conditions
# val1par: values of the changing parameter to evaluate
# nam1par: name of the changing parameter
# fmod: model
# tt: time point or range (tt=0 by default)
# ss: selection of ssi (0) or sse(1)
ssa1p <- function(initCon, parSet, mxInCon, val1par, nam1par, fmod, tt=0, ss){
  lsps <- lapply(1:dim(parSet)[1],function(i,parSet,val1par,nam1par,initCon,mxInCon,fmod,tt,ss){
    print("#############")
    print(i)
    print("#############")
    pp <- as.numeric(parSet[i,])
    names(pp) <- colnames(parSet)
    if(ss==0){
      ST1 <- lapply(val1par,function(p,nam1par,initCon,fmod){x<-pp;x[nam1par]<-p;ssi(x,nam1par,initCon,fmod,tt=0)},nam1par,initCon,fmod)
    }
    else if(ss==1){
      ST1 <- lapply(val1par,function(p,nam1par,initCon,fmod){x<-pp;x[nam1par]<-p;sse(x,nam1par,initCon,fmod,tt)},nam1par,initCon,fmod)
    }
    ST1df <-do.call(rbind,lapply(ST1,function(m){m <- m[which(!duplicated(round(m,10))),]}))
    ST1df <- ST1df[apply(ST1df[,- which(colnames(ST1df)%in%c(nam1par,"st"))],1, function(row){all(row<=max(mxInCon) & row>=min(mxInCon))}),]
    return(ST1df)
  }, parSet, val1par,nam1par,initCon,mxInCon,fmod,tt,ss)
  return (lsps)
}

#calculate bistable region
# ST1df: data frame containing the results of ssa1p
# mol: selected species to calculate the bistable region
# nam1par: name of the parameter that changes
# cutoffval: where to start to look for the min and max values of each steady state (by default 2)
BiOsreg <- function(ST1df, mol, nam1par, cutoffval=2){
  upbr <- ST1df[which(ST1df$st%in%c(2,4) & ST1df[,mol]>cutoffval & ST1df[,nam1par]>0),]
  dwbr <- ST1df[which(ST1df$st%in%c(2,4) & ST1df[,mol]<cutoffval & ST1df[,nam1par]>0),]
  print(c(dwbr[which.max(dwbr[,mol]), nam1par], upbr[which.min(upbr[,mol]), nam1par]))
  brp <- c(dwbr[which.max(dwbr[,mol]), nam1par], upbr[which.min(upbr[,mol]), nam1par])
  return(brp)
}

# remove misscalculations
# i: selected element (data frame) to remove the misscalculations
# lsps: the list of data frames that comes from ssa1p
# bsre: region of bistability or oscillations
# nam1par: name of the changing parameter
rmmiscal <- function(i,lsps,bsre,nam1par){
  mtx <- lsps[[i]][which(lsps[[i]][,nam1par]!=0),]
  if(length(bsre[[i]]) == 2){#if there is a bistable region
  ST1df <- do.call(rbind,lapply(unique(lsps[[i]][,nam1par]),function(x,m,r){#for unique values of namlpar, do:
    sbm <- m[m[,nam1par]==x,]#for all nam1par that are equal to the value we are evaluating
    if(all(x>=min(r) & x<=max(r))){# if value is in the bistable region
      if("3"%in%unique(sbm$st)){sbm3 <- unique(round(sbm[which(sbm$st==3),],2))}else{sbm3 <-NULL}
      if (any(c("2","4")%in%unique(sbm$st))){
        mx <- max(sbm[sbm$st%in%c(2,4),1])#max value for a molecule
        mn <- min(sbm[sbm$st%in%c(2,4),1])#min value for a molecule
        sbm0 <- apply(sbm[sbm$st==0,],1,function(row){if(all(row[1]<mx & row[1]>mn)){return(row)}})
        if (is.matrix(sbm0)){sbm0<-t(sbm0)}
        if (is.list(sbm0)){sbm0 <- do.call(rbind,sbm0)}
        sbm <- unique(round(rbind(sbm0,sbm[sbm$st%in%c(2,4),]),2))#round to 2 ss and select unique values of 2,4, 0
      }else{sbm <- NULL}
      sbm <- rbind(sbm, sbm3)
    }
    else if ("2"%in%unique(sbm$st)){
      sbm<-unique(round(sbm[which(sbm$st==2),],2))
      if (dim(sbm)[1]>1){
        #sbm <- rbind(sbm[which(sbm[,1]==min(sbm[,1])),], sbm[which(sbm[,1]==max(sbm[,1])),])
        sbm <- NULL
        }
    }
    else if("3"%in%unique(sbm$st)){sbm<-unique(round(sbm[which(sbm$st==3),],2))}
    else{sbm<-NULL}
    return(sbm)
  },mtx,bsre[[i]]))
  }else{ST1df <- mtx}#if there is not, return the same object
  #print((ST1df))
  return(ST1df)
}

# YT function to evaluate all time points and returns the max and min values
# x: time point to evaluate
# tc: output of ode evaluated for the specific model
YT <- function(x,tc){
  rx <- which(tc[,"time"]==x)#time point
  tp <- c(rx-1, rx, rx+1)#index for the time point (rx), the previous and the next one
  if(all(tp%in%seq(1,nrow(tc),1))){#all points in the range of time
    m <- tc[tp,-which(colnames(tc)=="time")]#values for the molecules
    if(any(is.na(m[,1]))){vm <- NULL}#if we don't have a value to compare, then NULL value is returned
    else{
      vm <- unlist(apply(m,2,function(cc){
        logval <- c(cc[1]<cc[2]&cc[2]>cc[3], cc[1]>cc[2]&cc[2]<cc[3]) #the evaluated value is max or min
        if( any(logval) ){return(cc[2])}else{}})) #the min or max is returned
    }
    if (!is.null(vm)){
      names(vm)<- gsub("[0-9]" , "", names(vm))
      names(vm)<- gsub(".", "", names(vm), fixed = T)
      vm <- data.frame(t=x,v=vm,n=names(vm),stringsAsFactors = F)#return time, value, molecule
      row.names(vm) <- NULL
      return(vm)
    }
  }#it is not in the range of time for searching
}

# ATfun calculates de amplitude and period for an oscillatory behaviour.
# parset: parameters
# timeseq: initial, final, step for time
# incon: initial conditions
# inframe: initial, final, step : seq(1.64,2.41,0.1)
# thrs: from this value, the step of the time sequence increases
# period: calculate the period: logical vector, T by default
ATfun <- function(inframe,parset,timeseq,incon,thrs,period=T,mod){ 
  AT <- lapply(inframe,function(p,pp,tt,ic,thrs,period,mod){#for each value of the selelected paramer (p)
    print(p)#print the parameter
    x <- pp; x[vn1]<-p #substitute it in the vector of parameters
    if (p <= thrs){tt <- c(0.0, 400,0.01)}#if the parameter is below the threshold, the step for time is smaller
    tc <-  ode(ic,seq(tt[1], tt[2], by = tt[3]), mod, x, maxsteps = 10000)#solve the model through ode function
    tc <- tc[which(tc[,"time"]>=(tt[2]/2)),]#take the second half for the rest of calculations (avoid the transient period)
    #calculate the amplitude through the min and max values
    Mnx <- t(do.call(rbind,lapply(colnames(tc)[-1],function(n,m){return(c(mn=min(m[,n]),mx=max(m[,n])))},tc)))
    colnames(Mnx) <- colnames(tc)[-1]
    #calculate the period. It also calculates the amplitude through the evaluation of all values
    if(period){
      #if there is no different between min and max, there is no period
      if(all(apply(Mnx,2,function(x){round(x[1],3)==round(x[2],3)}))){
        ytunq <- as.data.frame(do.call(rbind,lapply(colnames(Mnx),function(n,m){c(n=n,Ax=Mnx["mx",n],An=Mnx["mn",n],t=0)},Mnx)),stringsAsFactors = F)
      }
      #if there are differences, we calculate the period
      else{
        ytall <- sapply(tc[,"time"],YT,tc) #extract all min and max for each molecule. It returns time, value, molecule
        ytall <- do.call(rbind,ytall[-which(sapply(ytall, is.null))])#create a data frame with time (t), value (v), molecule (n)
        ytunq <- as.data.frame(do.call(rbind,lapply(unique(ytall[,"n"]),function(n,m){#extract min/max amplitue and period per molecule
          yt <- NULL
          m <- m[which(m$n==n),]#selecting values per molecule
          #calculate the Ax (max amplitude) and An (min amplitude)
          if(dim(m)[1] <= 1){yt <- NULL}
          else{yt <- c(n,max(round(as.numeric(m[,"v"]),3)),min(round(as.numeric(m[,"v"]),3)))}
          #calculate the period (t)
          if(period){
            if (dim(m)[1] <= 2){yt <- c(yt,0)}
            else{
              mp <- m[(nrow(m)-3):nrow(m),]
              if(round(as.numeric(mp[1,"v"],3))>round(as.numeric(mp[2,"v"],3))){yt <- c(yt,as.numeric(mp[3,"t"]) - as.numeric(mp[1,"t"]))}
              else if(round(as.numeric(mp[1,"v"],3))<round(as.numeric(mp[2,"v"],3))){yt <- c(yt,as.numeric(mp[4,"t"]) - as.numeric(mp[2,"t"]))}
              else{yt <- c(yt,0)}
            }
          }else{yt <- c(yt,NA)}
          
          return(yt)
        },ytall)), stringsAsFactors = F)
      }
      # if we don't want to calculate the period
    }else{
      ytunq <- as.data.frame(do.call(rbind,lapply(colnames(Mnx),function(n,m){c(n=n,Ax=Mnx["mx",n],An=Mnx["mn",n],t=0)},Mnx)),stringsAsFactors = F)
    }
    if(ncol(ytunq) != 4){}
    else{
      colnames(ytunq) <- c("n", "Ax", "An", "t")
      ytunq$p <- rep(x[vn1])
      return(ytunq)
    }
  },parset,timeseq,incon,thrs,period,mod)
  ATdf <- as.data.frame(do.call(rbind, AT),stringsAsFactors=F)
  ATdf$An <- as.numeric(ATdf$An)
  ATdf$Ax <- as.numeric(ATdf$Ax)
  ATdf[,"t"] <- as.numeric(ATdf[,"t"])
  
  return(ATdf)
}
