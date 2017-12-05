# Z matrix (internal use only)
Zmat <- function(x,type,theta) {
  #####
  lagFun <- function(x,maxlag) {
    if(maxlag>0) {
      out <- x
      if(maxlag>0) {
        for(w in 1:maxlag) {
          z <- rep(NA,length(x))
          if(w<length(x)) {
            for(i in (1+w):length(x)) {
              z[i] <- x[i-w]
              }
            }
          out <- cbind(out,z)        
          }
        }
      colnames(out) <- NULL
      out
      } else {
      rep(NA,length(x))
      }
    }
  #####
  if(type=="quec") {
    minlag <- theta[1]
    maxlag <- theta[2]
    laglim <- min(maxlag,length(x)-3)
    H <- c()
    for(i in 0:laglim) {
      if(i<minlag | i>maxlag) {
        H[i+1] <- 0
        } else {
        H[i+1] <- -4/(maxlag-minlag+2)^2*(i^2-(minlag+maxlag)*i+(minlag-1)*(maxlag+1))
        }
      }
    } else if(type=="qdec") {
    minlag <- theta[1]
    maxlag <- theta[2]
    laglim <- min(maxlag,length(x)-3)
    H <- c()
    for(i in 0:laglim) {
      if(i<minlag | i>maxlag) {
        H[i+1] <- 0
        } else {
        H[i+1] <- (i^2-2*maxlag*i+maxlag^2)/(maxlag-minlag)^2
        }
      }
    } else if(type=="gamma") {
    delta <- theta[1]
    lambda <- theta[2]
    xa <- 0  #####
    xb <- qgamma(0.99,1/(1-delta),-log(lambda))+xa
    laglim <- min(xb,length(x)-3)
    H <- c()
    for(i in 0:laglim) {
      if(i>=xa) {
        bnum <- (i-xa+1)^(delta/(1-delta))*lambda^(i-xa)
        xM <- (delta/(delta-1))/log(lambda)+xa-1
        bden <- (xM-xa+1)^(delta/(1-delta))*lambda^(xM-xa)
        H[i+1] <- bnum/bden
        } else {
        H[i+1] <- 0     
        }
      }                             
    }
  lagFun(x,laglim)%*%matrix(H,ncol=1)
  }

# check if a variable is quantitative (internal use only)
isQuant <- function(x) {
  if(is.numeric(x)) {
    T
    } else {
    if(!is.factor(x) && sum(!is.na(x))==0) {
      T
      } else {
      F
      }
    }
  }

# compute the lag limit (internal use only)
findLagLim <- function(data,group=NULL) {
  if(is.null(group)) {
    nrow(na.omit(data))
    } else {
    auxlaglim <- c()
    gruppi <- unique(data[,group])
    for(i in 1:length(gruppi)) {
      auxlaglim[i] <- nrow(na.omit(data[which(data[,group]==gruppi[i]),]))
      }
    min(auxlaglim,na.rm=T)
    }
  }

# adjust qualitative variables with less than two states (internal use only)
adjX <- function(X) {
  Xdel <- na.omit(X)
  for(i in 1:ncol(X)) {
    if(isQuant(X[,i])==F && length(unique(Xdel[,i]))<2) {    
      X[,i] <- rep(0,nrow(X))
      }
    }
  X
  }

# compute model score (internal use only)
modscore <- function(x,method) {
  res <- residuals(x)
  yfit <- fitted(x)
  n <- length(res)
  s2 <- summary(x)$sigma^2
  ll <- -0.5*n*log(2*pi*s2)-0.5/s2*sum(res^2)
  npar <- 2*(length(x$coefficients)+1)
  if(method=="ll") {
    ll
    } else if(method=="aic") {
    -2*ll+2*npar
    } else if(method=="bic") {
    -2*ll+log(n)
    } else if(method=="mdl") {
    (n-npar)*log(sum(res^2)/npar)+npar*log(sum(yfit^2))+(n-npar-1)*log(n/(n-npar))-(npar+1)*log(npar)
    } else {
    stop("Invalid method. Chose one among 'll', 'aic', 'bic' and 'mdl'",call.=F)
    }
  }

# least squares estimate (internal use only)
dlagLS <- function(y,X,group,data,type,theta,L) {
  if(is.null(group)) {
    xtime <- 1:nrow(data)
    Z <- data.frame(rep(NA,nrow(data)))
    auxcall <- c()
    Znam <- y
    if(length(X)>0) {
      for(i in 1:length(X)) {
        if(type[i]=="none") {
          iZ <- data[,X[i]]
          } else {
          iZ <- Zmat(data[,X[i]],type[i],theta[[i]])
          }
        Znam <- c(Znam,X[i])
        Z <- cbind(Z,iZ)
        }
      colnames(Z) <- Znam
      Z[,1] <- data[,y]
      form0 <- paste(y," ~ ",paste(colnames(Z)[-1],collapse="+"),sep="")
      } else {
      Z <- data.frame(data[,y])
      colnames(Z) <- y
      form0 <- paste(y," ~ 1",sep="")
      }
    Z <- adjX(Z)
    if(nrow(na.omit(Z))<3) stop("Insufficient sample size. Try to reduce lag leads",call.=F)
    mod0 <- lm(formula(form0),data=Z)
    } else {
    allZ <- xtime <- c()
    gruppi <- unique(data[,group])
    for(w in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[w])
      xtime[auxind] <- 1:length(auxind)
      Z <- data.frame(rep(NA,length(auxind)))
      auxcall <- c()
      Znam <- y
      if(length(X)>0) {
        for(i in 1:length(X)) {
          if(type[i]=="none") {
            iZ <- data[auxind,X[i]]
            } else {
            iZ <- Zmat(data[auxind,X[i]],type[i],theta[[i]])
            }
          Znam <- c(Znam,X[i])
          Z <- cbind(Z,iZ)
          }
        }
      colnames(Z) <- Znam
      Z[,1] <- data[auxind,y]
      allZ <- rbind(allZ,Z)
      }
    allZ <- cbind(data[,group],allZ)
    colnames(allZ)[1] <- group
    if(length(X)>0) {
      form0 <- paste(y," ~ -1+",group,"+",paste(colnames(allZ)[-c(1:2)],collapse="+"),sep="")
      } else {
      form0 <- paste(y," ~ -1+",group,sep="")
      }                               
    allZ <- adjX(allZ)
    if(nrow(na.omit(allZ))<3) stop("Insufficient sample size. Try to reduce lag leads",call.=F)
    mod0 <- lm(formula(form0),data=allZ)
    }
  bhat <- coefficients(mod0)
  bvar <- summary(mod0)$cov.unscaled
  auxdim <- ncol(bvar)
  bNA <- names(which(is.na(bhat)))
  if(length(bNA)>0) {
    bvar <- cbind(bvar,matrix(NA,nrow=nrow(bvar),ncol=length(bNA)))
    bvar <- rbind(bvar,matrix(NA,nrow=length(bNA),ncol=ncol(bvar)))
    colnames(bvar)[(auxdim+1):ncol(bvar)] <- rownames(bvar)[(auxdim+1):nrow(bvar)] <- bNA
    bvar <- bvar[names(bhat),names(bhat)]
    }
  #if(L>0) {
  #  res0 <- residuals(mod0)
  #  W <- matrix(nrow=length(res0),ncol=length(res0))
  #  for(i in 1:length(res0)) {
  #    for(j in 1:length(res0)) {
  #      ijlag <- abs(xtime[as.numeric(names(res0)[i])]-xtime[as.numeric(names(res0)[j])])
  #      W[i,j] <- max(0,(1-ijlag/(1+L)))*res0[i]*res0[j]
  #      }
  #    }
  #  Xmat <- model.matrix(mod0)
  #  Imat <- solve(t(Xmat)%*%Xmat)
  #  mod0$vcov <- Imat%*%t(Xmat)%*%W%*%Xmat%*%Imat
  #  } else {
  #  mod0$vcov <- bvar*summary(mod0)$sigma^2
  #  }
  #mod0$panel <- panelList
  mod0
  }

# scan formula (internal use only)
scanForm <- function(x) {
  auxform <- gsub(" ","",formula(x))[-1]
  ynam <- gsub(" ","",auxform[1])
  auX <- gsub(" ","",strsplit(auxform[2],"\\+")[[1]])
  if(sum(!is.na(auX)==0)) auX <- "1"
  auX <- setdiff(auX,"1")
  if(length(auX)>0) {
    lnames <- ltype <- rep(NA,length(auX))
    lpar <- list()
    for(i in 1:length(auX)) {        
      if(nchar(auX[i])>0) {
        if(identical("quec(",substr(auX[i],1,5))) {                           
          istr <- gsub("\\)","",gsub("quec\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "quec"
          lpar[[i]] <- as.numeric(istr[-1])
          } else if(identical("qdec(",substr(auX[i],1,5))) {                            
          istr <- gsub("\\)","",gsub("qdec\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "qdec"
          lpar[[i]] <- as.numeric(istr[-1])             
          } else if(identical("gamma(",substr(auX[i],1,6))) {
          istr <- gsub("\\)","",gsub("gamma\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "gamma"
          lpar[[i]] <- as.numeric(istr[-1])             
          } else {
          lnames[i] <- auX[i]
          ltype[i] <- "none"
          lpar[[i]] <- NA
          }                        
        if(grepl("\\(",lnames[i])) {
          stop("Unknown lag shape: ",strsplit(lnames[i],"\\(")[[1]][1],call.=F)
          }
        }    
      }
    names(lpar) <- names(ltype) <- lnames
    } else {
    lpar <- ltype <- c()
    }
  list(y=ynam,X=auX,ltype=ltype,lpar=lpar)
  }

# gamma constraints (internal use only)
gammaCons <- list()
gamseq <- seq(0.01,0.99,by=0.01)
gammaCons[[1]] <- gammaCons[[2]] <- matrix(nrow=length(gamseq),ncol=length(gamseq))
rownames(gammaCons[[1]]) <- colnames(gammaCons[[1]]) <- rownames(gammaCons[[2]]) <- colnames(gammaCons[[2]]) <- gamseq
names(gammaCons) <- c("sx","dx")
for(j in 1:length(gamseq)) {
  for(k in 1:length(gamseq)) {
    gammaCons[[1]][j,k] <- floor(qgamma(0.01,1/(1-gamseq[j]),-log(gamseq[k])))
    gammaCons[[2]][j,k] <- ceiling(qgamma(0.99,1/(1-gamseq[j]),-log(gamseq[k])))
    }
  }

# generate search grid (internal use only)
searchGrid <- function(maxgs,minwd,maxld,lag.type) {
  if(lag.type=="gamma") {
    auxind <- which(gammaCons$sx<=maxgs & gammaCons$dx<=maxld & gammaCons$dx-gammaCons$sx>=minwd,arr.ind=T)
    n <- nrow(auxind)
    if(n>0) {
      res <- c()
      for(i in 1:nrow(auxind)) {
        res <- rbind(res,as.numeric(c(rownames(gammaCons$sx)[auxind[i,1]],colnames(gammaCons$sx)[auxind[i,2]])))
        }
      outind <- seq(1,nrow(auxind),by=ceiling(n/500))
      res[outind,]
      }
    } else {
    auxmat <- c()
    for(i in 0:maxld) {
      for(j in i:maxld) { 
        if(i<=maxgs & j-i>=minwd)
        auxmat <- rbind(auxmat,c(i,j))
        }   
      }
    auxmat
    }
  }
  
# fit a distributed-lag linear regression model (internal use only)
dlaglm <- function(formula,group=NULL,data,L=0,adapt=FALSE,max.gestation=NULL,min.width=NULL,max.lead=NULL,sign=NULL,selection="aic") {
  auxscan <- scanForm(formula)
  y <- auxscan$y
  lagPar <- auxscan$lpar
  lagType <- auxscan$ltype
  lagNam <- names(lagPar)  
  if(length(lagNam)>0) {
    if(is.null(group)) {
      maxlag <- nrow(na.omit(data[,c(y,lagNam)]))-2
      } else {
      auxmlag <- c()
      gruppi <- unique(data[,group])
      for(i in 1:length(gruppi)) {
        auxind <- which(data[,group]==gruppi[i])
        if(length(auxind)>0) auxmlag[i] <- nrow(na.omit(data[auxind,c(y,lagNam)]))-2
        }
      maxlag <- min(auxmlag,na.rm=T)
      }
    } else {
    adapt <- F
    }
  laglimit <- findLagLim(data[,c(y,lagNam,group)],group=group)
  if(adapt==F) {
    modOK <- dlagLS(y=y,X=lagNam,group=group,data=data,type=lagType,theta=lagPar,L=L)
    if(length(c(group,auxscan$X))>0) {
      modOK$call$formula <- formula(paste(y," ~ ",paste(c(group,auxscan$X),collapse="+"),sep=""))
      } else {
      modOK$call$formula <- formula(paste(y," ~ 1",sep=""))
      }
    } else {
    bestPar <- vector("list",length=length(lagNam))
    names(bestPar) <- lagNam
    xOK <- c()
    fine <- 0
    while(fine==0) {
      xtest <- setdiff(lagNam,xOK)
      if(length(xtest)>0) {
        currentAIC <- rep(NA,length(xtest)) 
        currentP <- rep(NA,length(xtest))
        names(currentAIC) <- xtest
        names(currentP) <- xtest
        for(i in 1:length(xtest)) {
          if(xtest[i] %in% names(sign)) {
            isign <- sign[xtest[i]]
            } else {
            isign <- NULL
            }
          if(lagType[xtest[i]]!="none") {
            if(xtest[i] %in% names(max.gestation)) {
              maxgs <- max.gestation[xtest[i]]
              } else {
              maxgs <- laglimit
              }
            if(xtest[i] %in% names(min.width)) {
              minwd <- min.width[xtest[i]]
              } else {
              minwd <- 0
              }
            if(xtest[i] %in% names(max.lead)) {
              maxld <- max.lead[xtest[i]]
              } else {
              maxld <- laglimit
              }
            auxcons <- searchGrid(maxgs,minwd,maxld,lagType[xtest[i]])
            aic0 <- bhat0 <- c()
            pval0 <- c()
            for(j in 1:nrow(auxcons)) {
              testType <- lagType[c(xOK,xtest[i])]     
              testPar <- lagPar[c(xOK,xtest[i])]                       
              testPar[[xtest[i]]] <- auxcons[j,]
              mod0 <- dlagLS(y=y,X=names(testPar),group=group,data=data,type=testType,theta=testPar,L=0)
              summ0 <- summary(mod0)$coefficients
              ixstr <- paste(xtest[i],".",lagType[xtest[i]],sep="")                                                        
              if(ixstr %in% rownames(summ0)) {
                bhat0[j] <- summ0[ixstr,1]
                pval0[j] <- summ0[ixstr,4]  
                } else {
                bhat0[j] <- NA
                pval0[j] <- NA
                }                        
              aic0[j] <- modscore(mod0,selection)
              }        
            if(!is.null(isign)) {
              if(isign=="+") {
                auxsign <- which(bhat0>0)
                } else {
                auxsign <- which(bhat0<0)                  
                }
              if(length(auxsign)>0) {
                auxbest <- auxsign[which.min(aic0[auxsign])]
                } else {
                auxbest <- which.min(aic0)
                }
              } else {
              auxbest <- which.min(aic0)
              }
            bestPar[[xtest[i]]] <- auxcons[auxbest,]
            currentP[xtest[i]] <- pval0[auxbest]
            currentAIC[xtest[i]] <- aic0[auxbest]
            } else {
            testType <- lagType[c(xOK,xtest[i])]     
            testPar <- lagPar[c(xOK,xtest[i])]                       
            mod0 <- dlagLS(y=y,X=names(testPar),group=group,data=data,type=testType,theta=testPar,L=0)
            auxsumm <- summary(mod0)$coefficients
            if(xtest[i] %in% rownames(auxsumm)) {
              currentP[xtest[i]] <- auxsumm[xtest[i],4]          
              currentAIC[xtest[i]] <- modscore(mod0,selection)
              } else {
              currentP[xtest[i]] <- 1          
              currentAIC[xtest[i]] <- Inf            
              }
            }
          }
        issig <- which(currentP<0.05)
        if(length(issig)>0) {
          xOK <- c(xOK,names(currentAIC[issig])[which.min(currentAIC[issig])])
          } else {
          xOK <- c(xOK,names(currentAIC)[which.min(currentAIC)])          
          }
        } else {
        fine <- 1
        }
      }
    modOK <- dlagLS(y=y,X=lagNam,group=group,data=data,type=lagType,theta=bestPar,L=L)
    auxcall <- paste(y," ~ ",paste(c(group,lagNam),collapse="+"),sep="")
    for(i in 1:length(lagNam)) {
      if(lagType[i]!="none") auxcall <- gsub(lagNam[i],paste(lagType[i],"\\(",lagNam[i],",",paste(bestPar[[i]],collapse=","),"\\)",sep=""),auxcall)
      }
    modOK$call$formula <- formula(auxcall)
    }
  if(!is.null(group)) {
    grounam <- paste(group,sort(unique(data[,group])),sep="")
    auxcheck <- setdiff(grounam,names(modOK$coefficients)) 
    if(length(auxcheck)>0) {
      modOK$xlevels[[group]] <- gruppi
      modOK$coefficients <- c(rep(0,length(auxcheck)),modOK$coefficients)
      modOK$effects <- c(rep(0,length(auxcheck)),modOK$effects)
      names(modOK$coefficients)[1:length(auxcheck)] <- names(modOK$effects)[1:length(auxcheck)] <- auxcheck
      ##### modOK$qr
      }
    }
  modOK
  }

# compute lag effects of a covariate (internal use only)
lagEff <- function(model,x,cumul=F,conf=0.95,lag=NULL) {
  formstr <- strsplit(gsub(" ","",as.character(model$call$formula)[3]),"\\+")[[1]]
  auxscan <- scanForm(model$call)
  if(auxscan$ltype[x]=="quec") {
    sx <- auxscan$lpar[[x]][1]
    dx <- auxscan$lpar[[x]][2]
    imu <- model$coefficients[x]
    icov <- vcov(model)[x,x]
    if(is.null(lag)) {
      xgrid <- 0:dx
      } else {
      xgrid <- lag
      }
    iH <- c()
    for(i in 1:length(xgrid)) {
      if(xgrid[i]<sx-1 | xgrid[i]>dx+1) {
        iH[i] <- 0
        } else {
        iH[i] <- (-4/(sx-dx-2)^2)*(xgrid[i]^2-(sx+dx)*xgrid[i]+(sx-1)*(dx+1))
        }
      }
    iH <- matrix(iH,ncol=1)
    } else if(auxscan$ltype[x]=="qdec") {
    sx <- auxscan$lpar[[x]][1]
    dx <- auxscan$lpar[[x]][2]
    imu <- model$coefficients[x]
    icov <- vcov(model)[x,x]
    if(is.null(lag)) {
      xgrid <- 0:dx
      } else {
      xgrid <- lag
      }
    iH <- c()
    for(i in 1:length(xgrid)) {
      if(xgrid[i]<sx | xgrid[i]>dx) {
        iH[i] <- 0
        } else {
        iH[i] <- (xgrid[i]^2-2*(dx+1)*xgrid[i]+(dx+1)^2)/(dx-sx-2)^2
        }
      }
    iH <- matrix(iH,ncol=1)
    } else if(auxscan$ltype[x]=="gamma") {
    delta <- auxscan$lpar[[x]][1]
    lambda <- auxscan$lpar[[x]][2]
    sx <- 0  #####
    mingam <- qgamma(0.01,1/(1-delta),-log(lambda))+sx
    maxgam <- qgamma(0.99,1/(1-delta),-log(lambda))+sx
    imu <- model$coefficients[x]
    icov <- vcov(model)[x,x]
    if(is.null(lag)) {
      xgrid <- 0:maxgam
      } else {
      xgrid <- lag
      }
    iH <- c()
    for(i in 1:length(xgrid)) {
      if(xgrid[i]<mingam | xgrid[i]>maxgam) {
        iH[i] <- 0
        } else {
        bnum <- (xgrid[i]-sx+1)^(delta/(1-delta))*lambda^(xgrid[i]-sx)
        xM <- (delta/(delta-1))/log(lambda)+sx-1
        bden <- (xM-sx+1)^(delta/(1-delta))*lambda^(xM-sx)
        iH[i] <- bnum/bden
        }
      }
    iH <- matrix(iH,ncol=1)
    } else {  
    xgrid <- 0  
    imu <- model$coefficients[x]
    icov <- matrix(diag(rep(1,length(imu))),nrow=length(imu),ncol=length(imu))
    iH <- matrix(1,nrow=1,ncol=1)
    }                                    
  ibhat <- iH%*%imu
  ibse <- sqrt(diag(iH%*%icov%*%t(iH)))
  quan <- -qnorm((1-conf)/2)
  out <- cbind(ibhat,ibhat-quan*ibse,ibhat+quan*ibse)
  if(cumul==F) {
    rownames(out) <- xgrid
    colnames(out) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
    out
    } else {
    outC <- out
    for(j in 1:ncol(out)) {
      for(i in 1:nrow(out)) {
        outC[i,j] <- sum(out[1:i,j])
        }
      }
    rownames(outC) <- xgrid
    colnames(outC) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
    outC
    }
  }

# adf test (internal use only)
adft <- function(x,k=0) {
  k <- k + 1
  x <- as.vector(x,mode="double")
  y <- diff(x)
  n <- length(y)
  z <- embed(y,k)
  yt <- z[,1]
  xt1 <- x[k:n]
  tt <- k:n
  if(k>1) {
    yt1 <- z[,2:k]
    res <- lm(yt~xt1+1+tt+yt1)
    } else {
    res <- lm(yt~xt1+1+tt)
    }
  res.sum <- summary(res)
  STAT <- res.sum$coefficients[2,1]/res.sum$coefficients[2,2]
  table <- -1*cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), c(3.95, 
      3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 3.45, 3.43, 
      3.42, 3.41), c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), c(1.14, 
      1.19, 1.22, 1.23, 1.24, 1.25), c(0.8, 0.87, 0.9, 0.92, 
      0.93, 0.94), c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15, 
      0.24, 0.28, 0.31, 0.32, 0.33))
  tablen <- dim(table)[2]
  tableT <- c(25, 50, 100, 250, 500, 1e+05)
  tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
  tableipl <- numeric(tablen)
  for(i in (1:tablen)) {
    tableipl[i] <- approx(tableT,table[,i],n,rule=2)$y
    }
  PVAL <- approx(tableipl,tablep,STAT,rule=2)$y
  list(statistic=STAT,p.value=PVAL,'lag.order'=k-1)
  }

# kpss test (internal use only)
kpsst <- function(x,lshort=T) {
  x <- as.vector(x,mode="double")
  n <- length(x)
  t <- 1:n
  e <- residuals(lm(x~t))
  table <- c(0.216, 0.176, 0.146, 0.119)
  tablep <- c(0.01, 0.025, 0.05, 0.1)
  s <- cumsum(e)
  eta <- sum(s^2)/(n^2)
  s2 <- sum(e^2)/n
  if(lshort==T) {
    l <- trunc(3*sqrt(n)/13)
    } else {
    l <- trunc(10*sqrt(n)/14)
    }
  STAT <- eta/s2
  PVAL <- approx(table,tablep,STAT,rule=2)$y
  list(statistic=STAT,'trunc.lag'=l,p.value=PVAL)
  }

# apply differentiation (internal use only)
applyDiff <- function(x=NULL,group=NULL,time=NULL,data,k=0) {
  if(is.null(x)) x <- setdiff(colnames(data),c(group,time))
  deltaFun <- function(z,k) {
    if(k>0 & k<length(z)) {
      zd <- c(rep(NA,k),z[1:(length(z)-k)])
      z-zd
      } else if(k<=0) {
      z
      } else {
      rep(NA,length(z))
      }
    }
  diffdat <- data
  if(is.null(group)) {
    for(w in 1:length(x)) {
      if(isQuant(data[,x[w]])) {
        diffdat[,x[w]] <- deltaFun(data[,x[w]],k[w])      
        }
      }
    } else {
    data[,group] <- factor(data[,group])
    gruppi <- unique(data[,group])
    for(i in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[i])
      for(w in 1:length(x)) {
        if(isQuant(data[,x[w]])) {
          diffdat[auxind,x[w]] <- deltaFun(data[auxind,x[w]],k[w])
          }
        }
      } 
    }
  if(!is.null(time) & !is.null(group)) diffdat <- diffdat[order(data[,group]),]
  diffdat
  }

# unit root test
unirootTest <- function(x,group=NULL,time=NULL,data,test="adf",combine="choi",k=0,lshort=TRUE) {
  if(length(test)!=1 || (test %in% c("adf","kpss"))==F) stop("Argument 'test' must be either 'adf' or 'kpss'",call.=F)
  if(class(data)!="data.frame") stop("Argument 'data' must be a data.frame",call.=F)
  if(!is.null(group) && is.na(group)) group <- NULL
  if(!is.null(group) && length(group)!=1) stop("Argument 'group' must contain a single variable name",call.=F)
  auxvar <- setdiff(c(x,group,time),colnames(data))
  if(length(auxvar)>0) stop("Unknown variable: '",auxvar[1],"'",sep="",call.=F)
  if(length(x)!=1) stop("Argument 'x' must contain a single variable name",call.=F)
  if(isQuant(data[,x])==F) stop("'x' must be a numerical variable",call.=F)
  if(!is.null(time) && is.na(time)) time <- NULL
  if(!is.null(time) && length(time)!=1) stop("Argument 'time' must contain a single variable name",call.=F)
  if(!is.null(time) && isQuant(data[,time])==F) stop("'time' must be a numerical variable",call.=F)
  if(length(combine)!=1 || (combine %in% c("choi","demetrescu"))==F) stop("Argument 'combine' must be either 'choi' or 'demetrescu'",call.=F)
  if(length(k)!=1 || !is.numeric(k) || k<0 || k!=round(k)) stop("Argument 'k' must be an non-negative integer",call.=F)
  if(length(lshort)!=1 || !is.logical(lshort)) stop("Argument 'lshort' must be a logical value",call.=F)  
  ifelse(test=="adf", auxalt <- "stationary", auxalt <- "unit root")
  data[which(abs(data[,x])==Inf),x] <- NA
  options(warn=-1)
  if(is.null(group)) {                
    auxdat <- na.omit(data[,x])
    if(length(auxdat)>4 && var(auxdat)>0) {
      if(test=="adf") {
        auxadf <- adft(auxdat,k=k)
        } else {
        auxadf <- kpsst(auxdat,lshort=lshort)
        }
      res <- list(statistic=unname(auxadf$statistic),'lag.order'=k,alternative=auxalt,
        z.value=qnorm(auxadf$'p.value'),p.value=auxadf$'p.value',test=test,combine=NULL,n=length(auxdat))
      } else {
      res <- list(statistic=NULL,'lag.order'=k,alternative=auxalt,
        z.value=NULL,p.value=NULL,test=test,combine=NULL,n=length(auxdat))
      }
    } else {
    data[,group] <- factor(data[,group])
    gruppi <- unique(data[,group])
    auxstat <- auxp <- nwm <- c()
    for(i in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[i])
      auxdat <- data[auxind,x]
      nwm[i] <- sum(!is.na(auxdat))
      if(is.null(time)) {
        auxdat <- na.omit(auxdat)
        } else {
        auxdat <- na.omit(auxdat[order(data[,time])])        
        }                                                  
      if(length(auxdat)>4 && var(auxdat)>0) {
        if(test=="adf") {
          auxadf <- adft(auxdat,k=k)
          auxtest <- paste("Augmented Dickey-Fuller (lag order: ",k,")",sep="")
          } else {
          auxadf <- kpsst(auxdat,lshort=lshort)
          auxtest <- paste("KPSS (truncation parameter: ",auxadf$'trunc.lag',")",sep="")
          }
        auxstat[i] <- auxadf$statistic
        auxp[i] <- auxadf$p.value
        }
      }                       
    names(nwm) <- gruppi
    auxp <- na.omit(auxp)               
    auxstat <- na.omit(auxstat)  
    if(length(auxp)>0) {
      m <- length(auxp)
      rhat <- 1-var(qnorm(auxp))
      rstar <- max(rhat,-1/(m-1))
      if(combine=="demetrescu") {
        auxz <- sum(qnorm(auxp))/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
        } else if(combine=="choi") {
        auxz <- sum(qnorm(auxp))/sqrt(m)
        }
      auxpval <- 2*pnorm(-abs(auxz))
      res <- list(statistic=sum(nwm*auxstat)/sum(nwm),alternative=auxalt, #'lag.order'=k,
        z.value=auxz,p.value=auxpval,test=test,combine=combine,n=nwm)
      } else {
      res <- list(statistic=NULL,alternative=auxalt, #'lag.order'=k,
        z.value=NULL,p.value=NULL,test=test,combine=combine,n=nwm)
      }
    }
  options(warn=0)
  res
  }

# fit the imputation model (internal use only)
impuFit <- function(data,group=NULL,isLat) {
  xall <- setdiff(colnames(data),group)
  if(!is.null(group)) data[,group] <- factor(data[,group])
  xfact <- c()
  for(j in 1:length(xall)) {
    if(isQuant(data[,xall[j]])==F) xfact <- c(xfact,xall[j])
    }
  xcont <- setdiff(xall,xfact)
  res <- vector("list",length=length(xall))
  names(res) <- xall
  for(i in 1:length(xall)) {
    ipar <- c(group,setdiff(xall,xall[i]))
    if(length(ipar)>0) {
      iparstr <- paste(ipar,collapse="+")
      } else {
      iparstr <- "1"
      }
    if(xall[i] %in% xcont) {
      res[[xall[i]]] <- lm(formula(paste(xall[i],"~",iparstr,sep="")),data=data)  
      if(isLat[xall[i]]==1) {
        res[[xall[i]]]$residuals <- res[[xall[i]]]$residuals/sd(res[[xall[i]]]$residuals)
        }
      } else {
      res[[xall[i]]] <- multinom(formula(paste(xall[i],"~",iparstr,sep="")),data=data,trace=F)       
      #res[[xall[i]]] <- mlogit(formula(paste(xall[i],"~",iparstr,sep="")),data=data)       
      }
    }
  logL <- sum(unlist(lapply(res,logLik)))
  list(estimate=res,logL=logL)
  }

# missing imputation using EM (internal use only)
EM.imputation <- function(x,group=NULL,data,tol=0.0001,maxiter=500,quiet=FALSE) {
  xlev <- vector("list",length=ncol(data))
  for(i in 1:ncol(data)) {
    if(isQuant(data[,i])==F) {     
      data[,i] <- factor(data[,i])
      xlev[[i]] <- levels(data[,i])
      }
    }
  names(xlev) <- colnames(data)
  isLat <- rep(0,ncol(data))
  isNA <- matrix(0,nrow=nrow(data),ncol=ncol(data))
  colnames(isNA) <- names(isLat) <- colnames(data)
  n <- nrow(data)
  auxdat <- data
  for(i in 1:length(x)) {
    ina <- which(is.na(auxdat[,x[i]]))
    if(length(ina)>0){
      if(isQuant(auxdat[,x[i]])) {
        imu <- mean(data[,x[i]],na.rm=T)
        if(!is.na(imu)) {
          auxdat[ina,x[i]] <- rep(imu,length(ina))
          } else {
          auxdat[ina,x[i]] <- rnorm(length(ina))          
          isLat[x[i]] <- 1
          }
        } else {
        iauxM <- table(data[,x[i]])
        if(length(iauxM)<2) stop("Variable '",x[i],"' has less than two unique values and cannot be imputed",sep="",call.=F)
        imu <- names(iauxM)[which.max(iauxM)]
        auxdat[ina,x[i]] <- rep(imu,length(ina))               
        }
      isNA[ina,x[i]] <- 1
      }
    }
  fit0 <- impuFit(auxdat,group=group,isLat)
  currentMod <- fit0$estimate
  currentLik <- fit0$logL
  fine <- forcend <- 0
  count <- 1
  if(quiet==F) {
    cat("Starting EM...")
    flush.console() 
    }
  while(fine==0) {
    if(quiet==F) {
      cat('\r',"EM iteration ",count,". Log-likelihood: ",currentLik,sep="")
      flush.console() 
      }
    for(i in 1:length(x)) {
      if(sum(isNA[,x[i]])>0) {
        if(is.null(xlev[[x[i]]])) {
          auxdat[which(isNA[,x[i]]==1),x[i]] <- predict(currentMod[[x[i]]])[which(isNA[,x[i]]==1)]
          } else {
          ipred <- predict(currentMod[[x[i]]])[which(isNA[,x[i]]==1)]
          auxdat[which(isNA[,x[i]]==1),x[i]] <- xlev[[x[i]]][ipred]
          }
        }
      }
    newFit <- impuFit(auxdat,group=group,isLat)
    newMod <- newFit$estimate
    newLik <- newFit$logL
    if(newLik-currentLik>tol & count<maxiter) {
      currentMod <- newMod
      currentLik <- newLik
      count <- count+1
      } else {
      fine <- 1
      if(count>=maxiter) forcend <- 1
      }
    }
  if(quiet==F) {
    if(forcend==0) {
      cat('\r',"EM converged after ",count," iterations. Log-likelihood: ",newLik,sep="","\n")
      } else {
      cat('\r',"EM stopped after ",maxiter," iterations. Log-likelihood: ",newLik,sep="","\n")      
      }
    }
  auxdat
  }
 
# function to plot a lag shape (internal use only)
makeShape <- function(bmat,maxlag,cumul,conf,ylim,title) {
  if(!is.null(maxlag)) {
    ymlag <- max(as.numeric(rownames(bmat)))
    if(maxlag>=ymlag) {
      if(cumul==F) {
        addmat <- matrix(0,nrow=maxlag-ymlag+1,ncol=3)
        } else {
        addmat <- matrix(bmat[nrow(bmat),],nrow=maxlag-ymlag+1,ncol=3,byrow=T)          
        }
      bmat <- rbind(bmat,addmat)
      rownames(bmat) <- -1:(maxlag+1)
      } else {                
      bmat <- bmat[1:(which(rownames(bmat)==as.character(maxlag))+1),]
      }
    auxNZ <- which(bmat[,1]!=0)
    } else {
    auxNZ <- which(bmat[,1]!=0)
    if(cumul==F) {
      if(nrow(bmat)>max(auxNZ)+1) {
        bmat <- bmat[1:(max(auxNZ)+1),]
        }
      } else {                          
      auxCm <- min(intersect(which(diff(bmat[,1])==0)+1,auxNZ))
      if(nrow(bmat)>auxCm) {
        bmat <- bmat[1:auxCm,]
        }
      auxNZ <- which(bmat[,1]!=0)
      }
    }
  xaux <- (1:nrow(bmat))-2
  upLim <- 1.05*max(bmat)
  if(is.null(ylim)) {
    lowLim <- 1.05*min(bmat)
    upLim <- max(abs(c(upLim,lowLim)))
    lowLim <- -max(abs(c(upLim,lowLim)))
    } else {
    lowLim <- ylim[1]
    upLim <- ylim[2]
    }
  plot(0,type="n",xlim=c(min(xaux),max(xaux)),ylim=c(lowLim,upLim),yaxs="i",xaxs="i",cex.lab=1.2,
    lwd=2,xaxt="n",yaxt="n",xlab="Lag",ylab="Coefficient",main=title,cex.main=1.2) 
  if(cumul==T) {
    mtext("cumulative lag shape",cex=0.9)
    } else {
    mtext("instantaneous lag shape",cex=0.9)    
    }
  #####
  #yaxaux <- seq(lowLim,upLim,length=21)
  #ylabaux <- signif(yaxaux,3)
  #ylabaux[11] <- 0
  #xaxaux <- seq(min(xaux),max(xaux))
  #auxby <- max(1,round((max(xaux)-min(xaux)+1)/30))
  #xlabaux1 <- xlabaux2 <- seq(min(xaux),max(xaux),by=auxby)
  #xlabaux2[c(1,length(xlabaux1))] <- NA                                  
  #for(i in 1:length(auxNZ)) {
  #  ix <- xaux[auxNZ[i]]+c(-0.3,0.3)            
  #  iysx <- rep(bmat[auxNZ[i],2],2)      
  #  iydx <- rep(bmat[auxNZ[i],3],2)
  #  polygon(c(ix,rev(ix)),c(iysx,rev(iydx)),border=NA,col="grey80")
  #  }
  ##segments(xaux[auxNZ],bmat[auxNZ,2],xaux[auxNZ],bmat[auxNZ,3],col="grey35")
  ##segments(xaux[auxNZ]-0.1,bmat[auxNZ,2],xaux[auxNZ]+0.1,bmat[auxNZ,2],col="grey35")
  ##segments(xaux[auxNZ]-0.1,bmat[auxNZ,3],xaux[auxNZ]+0.1,bmat[auxNZ,3],col="grey35")
  #abline(h=yaxaux,v=seq(min(xaux),max(xaux),by=auxby),col="grey75",lty=2)
  #abline(h=0,lty=2,col="grey35")
  #lines(bmat[,1]~xaux,col="grey35",lty=2)
  #auxpoi <- max(1,min(auxNZ)-1):min(length(xaux),max(auxNZ)+1)
  #points(bmat[auxpoi,1]~xaux[auxpoi],col="grey35",lty=2,cex=0.6)
  #axis(1,at=xlabaux1,labels=xlabaux2,cex.axis=1.1)
  #axis(2,at=yaxaux,labels=ylabaux,cex.axis=1.1) 
  #####
  usesplin <- 0
  auxsplin <- bmat[auxNZ,1]
  if(length(auxsplin)>=3) usesplin <- 1 
  isTrunc <- 0
  auxM <- which.max(abs(bmat[,1]))
  lM <- as.numeric(rownames(bmat)[auxM])  
  if(auxM==1 || sum(bmat[1:(auxM-1),1])==0) isTrunc <- 1                                                                                                               
  if(usesplin==0) {
    xgrid <- xaux
    ygrid <- bmat
    } else {
    if(isTrunc==1) {
      xgrid <- seq(lM,max(xaux[auxNZ])+1,length=1000)
      splind <- auxM:min(length(xaux),(max(auxNZ)+1)) 
      } else {
      xgrid <- seq(min(xaux[auxNZ])-1,max(xaux[auxNZ])+1,length=1000)
      splind <- max(1,(min(auxNZ)-1)):min(length(xaux),(max(auxNZ)+1))
      }
    ygrid <- matrix(nrow=length(xgrid),ncol=ncol(bmat))
    for(i in 1:ncol(bmat)) {
      auxspl <- interpSpline(xaux[splind],bmat[splind,i])
      ygrid[,i] <- predict(auxspl,xgrid)$y
      }
    }
  if(isTrunc==1) {
    xgrid <- c(lM-1,xgrid)
    ygrid <- rbind(rep(0,3),ygrid)
    }
  #for(i in 1:length(auxNZ)) {
  #  ix <- xaux[auxNZ[i]]+c(-0.25,0.25)            
  #  iysx <- rep(bmat[auxNZ[i],2],2)      
  #  iydx <- rep(bmat[auxNZ[i],3],2)
  #  polygon(c(ix,rev(ix)),c(iysx,rev(iydx)),border=NA,col="grey80")
  #  }
  polygon(c(xgrid,rev(xgrid)),c(ygrid[,2],rev(ygrid[,3])),border=NA,col="grey80")
  yaxaux <- seq(lowLim,upLim,length=21)
  ylabaux <- signif(yaxaux,3)
  ylabaux[11] <- 0
  xaxaux <- seq(min(xaux),max(xaux))
  auxby <- max(1,round((max(xaux)-min(xaux)+1)/30))
  xlabaux1 <- xlabaux2 <- seq(min(xaux),max(xaux),by=auxby)
  xlabaux2[c(1,length(xlabaux1))] <- NA
  abline(h=yaxaux,v=seq(min(xaux),max(xaux),by=auxby),col="grey75",lty=2)
  abline(h=0,lty=2,col="grey35")                                        
  lines(ygrid[,1]~xgrid,col="grey40",lty=2)
  auxpoi <- max(1,min(auxNZ)-1):min(length(xaux),max(auxNZ)+1) ###
  points(bmat[auxpoi,1]~xaux[auxpoi],col="grey35",lty=2,cex=0.6)
  axis(1,at=xlabaux1,labels=xlabaux2,cex.axis=1.1)
  axis(2,at=yaxaux,labels=ylabaux,cex.axis=1.1)  
  #####
  if(cumul==F) {
    tcaVal <- signif(apply(bmat,2,sum),5)
    } else {
    tcaVal <- signif(bmat[nrow(bmat),],5)
    }
  confLeg <- paste("   ",conf*100,"% CI: (",tcaVal[2],", ",tcaVal[3],")",sep="")      
  if(max(bmat[,1])>0) {
    legpos <- "bottomright"
    } else {
    legpos <- "topright"
    }                                       
  est <- bmat[,1]
  if(cumul==T) {
    newest <- est[1]
    for(i in 2:length(est)) {
      newest[i] <- est[i]-est[i-1]
      } 
    est <- newest
    }                                       
  minlag <- min(as.numeric(rownames(bmat)[which(est!=0)]))
  maxlag <- max(as.numeric(rownames(bmat)[which(est!=0)]))                                                                             
  legend(legpos,legend=c(paste("Effective lags: ",minlag," to ",maxlag,sep=""),paste("Cumulative coefficient: ",tcaVal[1],sep=""),confLeg),bty="n",cex=1.1)
  box()
  }    

# plot the lag shape associated to an overall causal effect or a path
lagPlot <- function(x,from=NULL,to=NULL,path=NULL,maxlag=NULL,cumul=FALSE,conf=0.95,use.ns=FALSE,ylim=NULL,title=NULL) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(!is.null(maxlag) && (length(maxlag)!=1 || !is.numeric(maxlag) || maxlag<=0 || maxlag!=round(maxlag))) stop("Argument 'maxlag' must be a positive integer number",call.=F)
  if(length(cumul)!=1 || !is.logical(cumul)) stop("Arguent 'cumul' must be a logical value",call.=F)
  if(length(use.ns)!=1 || !is.logical(use.ns)) stop("Argument 'use.ns' must be a logical value",call.=F)
  if(!is.null(ylim) && (length(ylim)!=2 || ylim[1]>=ylim[2])) stop("Invalid argument 'ylim'",call.=F)
  if(!is.null(from) && !is.na(from) && (is.null(to) & is.null(path))) {
    path <- from
    auxstr <- strsplit(path,"\\*")[[1]]
    if(length(auxstr)<2) {
      stop("Argument 'to' is missing",call.=F)
      } else {  
      from <- NULL
      }
    }
  if(is.null(path) && (is.null(from) || is.na(from))) stop("Argument 'from' is missing",call.=F)
  if(is.null(path) && (is.null(to) || is.na(to))) stop("Argument 'to' is missing",call.=F)
  xedgF <- edgeMat(x,conf=conf,full=T)
  xedg <- edgeMat(x,conf=conf,full=F)
  if(is.null(path)) {                        
    auxpa <- causalEff(x,from=from,to=to,lag=NULL,cumul=cumul,conf=conf,use.ns=use.ns)$overall
    } else {
    auxstr <- strsplit(path,"\\*")[[1]]
    from <- auxstr[1]
    to <- rev(auxstr)[1] 
    isIn <- isInF <- 1
    for(i in 2:length(auxstr)) {
      isIn <- isIn*length(which(xedg[,1]==auxstr[i-1] & xedg[,2]==auxstr[i]))
      isInF <- isInF*length(which(xedgF[,1]==auxstr[i-1] & xedgF[,2]==auxstr[i]))
      } 
    if((use.ns==T & isInF>0) | isIn>0) {
      auxpa <- causalEff(x,from=from,to=to,cumul=cumul,conf=conf,use.ns=use.ns)     
      auxpa <- auxpa[[path]] 
      } else {
      if(isInF>0) {
        stop("Path not found. Try to reduce 'conf' or to set 'use.ns' to TRUE",call.=F)
        } else {
        stop("Inexistent path",call.=F)
        }
      }
    }               
  yaux <- rbind(rep(0,3),auxpa)
  rownames(yaux) <- c(-1:(nrow(yaux)-2))
  if(is.null(title)) {
    if(is.null(path)) {
      title <- paste(to," ~ ",from,sep="")
      } else {
      title <- paste(auxstr,collapse=" * ")
      }
    }
  makeShape(yaux,maxlag=maxlag,cumul=cumul,conf=conf,ylim=ylim,title=title)
  }

# check if a vector is named (internal use only)
hasname <- function(x) {
  out <- sum(nchar(names(x))>0)==length(x)
  if(is.na(out)) out <- F
  out
  }
  
# check control options (internal use only)
controlCheck <- function(parstr,control,parSets,is.sign=F) {
  if(parstr %in% names(control)) {
    if(!is.list(control[[parstr]]) || hasname(control[[parstr]])==F) {
      stop("Component '",parstr,"' of argument 'control' is not a named list",call.=F)
      }
    for(i in 1:length(control[[parstr]])) {
      if(hasname(control[[parstr]][[i]])==F) stop("Component '",parstr,"' of argument 'control' is not a named list: ",names(control[[parstr]])[i],call.=F)
      }
    auxch <- setdiff(names(control[[parstr]]),names(parSets))
    if(length(auxch)>0) stop("Unknown variable in component '",parstr,"' of argument 'control': ",auxch[1],call.=F)
    auxG <- lapply(control[[parstr]],names)
    gnam <- names(auxG)
    for(i in 1:length(gnam)) {
      auxpar <- parSets[[gnam[i]]]
      auxev <- auxG[[gnam[i]]] 
      for(j in 1:length(auxev)) {
        if((auxev[j] %in% auxpar)==F) stop("Invalid covariate in component '",parstr,"' of argument 'control': ",gnam[i]," | ",auxev[j],call.=F)
        }
      }
    for(i in 1:length(control[[parstr]])) { 
      ival <- control[[parstr]][[i]]
      for(j in 1:length(ival)) {
        if(is.sign==F) {
          if(!is.numeric(ival[j]) || ival[j]<0 || ival[j]!=round(ival[j])) {
            stop("Invalid control options in component '",parstr,"' of argument 'control': ",names(control[[parstr]])[i]," | ",names(control[[parstr]][[i]])[j],call.=F)
            }
          } else {
          if((ival[j] %in% c("+","-"))==F) {
            stop("Invalid control options in component '",parstr,"' of argument 'control': ",names(control[[parstr]])[i]," | ",names(control[[parstr]][[i]])[j],call.=F)
            }
          }
        }
      }
    }
  }

# fit a distributed-lag linear SEM
dlsem <- function(model.code,group=NULL,exogenous=NULL,data,log=FALSE,diff.options=list(test="adf",combine="choi",k=0,lshort=TRUE,maxdiff=3),
  imput.options=list(tol=0.0001,maxiter=500,no.imput=NULL),global.control=NULL,local.control=NULL) {
  #
  if(!is.list(model.code) || sum(sapply(model.code,class)!="formula")>0) stop("Argument 'model code' must be a list of formulas",call.=F)
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be an object of class 'data.frame'",call.=F)  
  nameOfData <- deparse(substitute(data))
  if(!is.null(group) && is.na(group)) group <- NULL
  if(!is.null(exogenous) && is.na(exogenous)) exogenous <- NULL
  if(!is.null(group) && length(group)!=1) stop("Argument 'group' must contain a single variable name",call.=F)
  if(!is.null(group) && (group %in% colnames(data))==F) stop("Variable '",group,"' not found in data",sep="",call.=F)
  if(length(log)!=1 || !is.logical(log)) stop("Argument 'log' must be a logical value",call.=F)
  if(!is.null(diff.options) && !is.list(diff.options)) stop("Argument 'diff.options' must be a list",call.=F)
  if(!is.null(imput.options) && !is.list(imput.options)) stop("Argument 'imput.options' must be a list",call.=F)
  if(!is.null(global.control) && !is.list(global.control)) stop("Argument 'global.control' must be a list",call.=F) 
  if(!is.null(local.control) && !is.list(local.control)) stop("Argument 'local.control' must be a list",call.=F)  
  #
  if(is.null(diff.options)) diff.options <- list(test="adf",combine="choi",k=0,lshort=T,maxdiff=3)
  auxch <- setdiff(names(diff.options),c("test","combine","k","lshort","maxdiff"))
  if(length(auxch)>0) stop("Invalid component '",auxch[1],"' in argument 'diff.options'",sep="",call.=F)
  if("test" %in% names(diff.options)) {
    test <- diff.options$test
    } else {
    test <- "adf"
    }
  if("combine" %in% names(diff.options)) {
    combine <- diff.options$combine
    } else {
    combine <- "choi"
    }
  if("k" %in% names(diff.options)) {
    k <- diff.options$k
    } else {
    k <- 0
    }
  if("lshort" %in% names(diff.options)) {
    lshort <- diff.options$lshort
    } else {
    lshort <- T
    }
  if("maxdiff" %in% names(diff.options)) {
    maxdiff <- diff.options$maxdiff
    } else {
    maxdiff <- 3
    }  
  if(length(test)!=1 || (test %in% c("adf","kpss"))==F) stop("Argument 'test' must be either 'adf' or 'kpss'",call.=F)
  if(length(combine)!=1 || (combine %in% c("choi","demetrescu"))==F) stop("Argument 'combine' must be either 'choi' or 'demetrescu'",call.=F)
  if(length(maxdiff)!=1 || !is.numeric(maxdiff) || maxdiff<0 || maxdiff!=round(maxdiff)) stop("Argument 'maxdiff' must be a non-negative integer number",call.=F)
  if(length(k)!=1 || !is.numeric(k) || k<0 || k!=round(k)) stop("Argument 'k' must be an non-negative integer",call.=F)
  if(length(lshort)!=1 || !is.logical(lshort)) stop("Argument 'lshort' must be a logical value",call.=F)  
  #
  if(is.null(imput.options)) imput.options <- list(tol=0.0001,maxiter=500,no.imput=NULL)
  auxch <- setdiff(names(imput.options),c("tol","maxiter","no.imput"))
  if(length(auxch)>0) stop("Invalid component '",auxch[1],"' in argument 'imput.options'",sep="",call.=F)
  if("tol" %in% names(imput.options)) {
    tol <- imput.options$tol
    } else {
    tol <- 0.0001
    }
  if("maxiter" %in% names(imput.options)) {
    maxiter <- imput.options$maxiter
    } else {
    maxiter <- 500
    }
  if("no.imput" %in% names(imput.options)) {
    no.imput <- imput.options$no.imput
    } else {
    no.imput <- NULL
    }
  if(length(maxiter)!=1 || maxiter<0 || maxiter!=round(maxiter)) stop("Argument 'maxiter' must be a non-negative integer number",call.=F)
  if(length(tol)!=1 || tol<=0) stop("Argument 'tol' must be a positive real number",call.=F)              
  #
  if(is.null(global.control)) {
    global.control <- list(adapt=F,selection="aic")
    } else {
    if(("adapt" %in% names(global.control))==F) global.control$adapt <- F
    if(("selection" %in% names(global.control))==F) global.control$selection <- "aic"
    }                                              
  auxch <- setdiff(names(global.control),c("adapt","max.gestation","min.width","max.lead","sign","selection"))
  if(length(auxch)>0) stop("Invalid component '",auxch[1],"' in argument 'global.options'",sep="",call.=F)  
  if(!is.logical(global.control$adapt)) stop("Component 'adapt' in argument 'global.options' must be a logical value",call.=F)
  if((global.control$selection %in% c("aic","bic","mdl"))==F) stop("Component 'selection' in argument 'global.options' must be one among 'aic', 'bic' and 'mdl'",call.=F)
  if(!is.null(global.control$max.gestation) && (global.control$max.gestation<0 || global.control$max.gestation!=round(global.control$max.gestation))) {
    stop("Component 'max.gestation' in argument 'global.options' must be a non-negative integer value",call.=F)
    }
  if(!is.null(global.control$min.width) && (global.control$min.width<0 || global.control$min.width!=round(global.control$min.width))) {
    stop("Component 'min.width' in argument 'global.options' must be a non-negative integer value",call.=F)
    }
  if(!is.null(global.control$max.lead) && (global.control$max.lead<0 || global.control$max.lead!=round(global.control$max.lead))) {
    stop("Component 'max.lead' in argument 'global.options' must be a non-negative integer value",call.=F)
    }
  if(!is.null(global.control$sign) && ((global.control$sign %in% c("+","-"))==F)) {
    stop("Component 'sign' in argument 'global.options' must be either '+' or '-'",call.=F)
    }
  if(!is.null(global.control$min.width) && !is.null(global.control$max.lead) && global.control$min.width>global.control$max.lead) {
    stop("Component 'min.width' greater than component 'max.lead' in argument 'global.control'",call.=F)
    }
  #
  rownames(data) <- 1:nrow(data)
  res <- pset <- list()
  messlen <- 0
  for(i in 1:length(model.code)) {
    icheck <- scanForm(model.code[[i]])
    ipar <- names(icheck$lpar)
    if(length(ipar)>0) {
      if(!is.null(group)) {
        if(group %in% ipar) stop("Variable '",group,"' is defined as a group factor and appears in 'model.code'",call.=F) 
        }
      auxexo <- intersect(ipar,exogenous)
      if(length(auxexo)>0) stop("Variable '",auxexo[1],"' appears both in 'model.code' and in 'exogenous'",call.=F)
      auxdupl <- duplicated(ipar)
      if(sum(auxdupl)>0) stop("Duplicated covariate for ",icheck$y,": ",ipar[auxdupl][1],call.=F)
      for(j in 1:length(icheck$lpar)) {
        if(icheck$ltype[j]=="quec") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>icheck$lpar[[j]][2]) {
            stop("Invalid settings for lag shape 'quec' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="qdec") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>icheck$lpar[[j]][2]) {
            stop("Invalid settings for lag shape 'qdec' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="gamma") {
          if(length(icheck$lpar[[j]])!=2 || sum(icheck$lpar[[j]]<=0)>0 || sum(icheck$lpar[[j]]>=1)>0) {
            stop("Invalid settings for lag shape 'gamma' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          }
        }
      }
    if(sum(grepl("\\-",model.code[[i]]))>0) stop("Invalid character '-' in 'model.code'",call.=F)
    if(sum(grepl("\\:",model.code[[i]]))>0) stop("Invalid character ':' in 'model.code'",call.=F) #####
    if(sum(grepl("\\*",model.code[[i]]))>0) stop("Invalid character '*' in 'model.code'",call.=F) #####
    auxstr <- as.character(model.code[[i]])[-1]
    ynam <- gsub(" ","",auxstr[1])
    names(model.code)[i] <- ynam
    xnam <- gsub(" ","",strsplit(auxstr[2],"\\+")[[1]])  
    xnam <- xnam[which(nchar(xnam)>0)]
    for(j in 1:length(xnam)) {
      if(grepl("\\(",xnam[j])) xnam[j] <- strsplit(strsplit(xnam[j],",")[[1]][1],"\\(")[[1]][2]
      }
    if(identical(xnam,"1") | length(xnam)==0) {
      pset[[i]] <- character(0)
      } else {
      pset[[i]] <- xnam
      }
    names(pset)[i] <- ynam 
    }   
  codnam <- c()
  for(i in 1:length(model.code)) {
    codnam[i] <- as.character(model.code[[i]])[2]
    }
  auxdupl <- duplicated(codnam)
  if(sum(auxdupl)>0) stop("Duplicated response variable: ",codnam[auxdupl][1],call.=F)
  nodenam <- unique(c(names(pset),unlist(pset)))
  auxadd <- setdiff(nodenam,codnam)
  if(length(auxadd)>0) {
    for(i in length(auxadd):1) {
      model.code <- c(formula(paste(auxadd[i],"~1",sep="")),model.code)
      names(model.code)[1] <- auxadd[i]
      pset[[length(pset)+1]] <- character(0) 
      names(pset)[length(pset)] <- auxadd[i]
      }
    }
  if(length(nodenam)<2) stop("The model cannot contain less than 2 endogenous variables",call.=F)
  auxvar <- setdiff(c(nodenam,exogenous),colnames(data))
  if(length(auxvar)>0) {
    stop("Variable '",auxvar[1],"' not found in data",sep="",call.=F)
    #data[,auxvar] <- as.numeric(rep(NA,nrow(data)))
    #cat(paste("Unobserved variables: ",paste(auxvar,collapse=", "),sep=""),"\n")
    }
  for(i in 1:length(nodenam)) {
    if(isQuant(data[,nodenam[i]])==F) stop("Variable '",nodenam[i],"' is qualitative and cannot appear in 'model.code'",sep="",call.=F)
    if(length(na.omit(data[,nodenam[i]]))==0) stop("Variable '",nodenam[i],"' has no observed values in data",sep="",call.=F)
    }
  if(!is.null(exogenous)) {
    for(i in 1:length(exogenous)) {
      if(length(na.omit(data[,exogenous[i]]))==0) stop("Variable '",exogenous[i],"' has no observed values in data",sep="",call.=F)
      }                                                                                                                                                                          
    }
  #for(i in 1:length(model.code)) {
  #  icheck <- scanForm(model.code[[i]])
  #  ilimit <- findLagLim(data[,c(icheck$y,names(icheck$lpar),group)],group=group)
  #  }
  G <- new("graphNEL",nodes=names(pset),edgemode="directed")    
  for(i in 1:length(pset)) {
    if(length(pset[[i]])>0) {
      for(j in 1:length(pset[[i]])) {    
        G <- addEdge(pset[[i]][j],names(pset)[i],G,1) 
        }
      }
    }                   
  topG <- topOrder(G)
  if(is.null(topG)) stop("The model code implies directed cycles",call.=F)  
  auxch <- setdiff(no.imput,topG)
  if(length(auxch)>0) stop("Unknown variable '",auxch[1],"' in component 'no.imput' of argument 'imput.options'",call.=F)  
  auxch <- setdiff(names(local.control),c("adapt","max.gestation","min.width","max.lead","sign"))  #,"L"
  if(length(auxch)>0) stop("Unknown component '",auxch[1],"' in argument 'local.control'",call.=F)
  #if("L" %in% names(local.control)) {
  #  if(!is.numeric(local.control[["L"]]) || hasname(local.control[["L"]])==F) stop("Component 'L' of argument 'local.control' must be a named numeric vector",call.=F)   
  #  auxch <- setdiff(names(local.control[["L"]]),topG)
  #  if(length(auxch)>0) stop("Unknown variables in component 'L' of argument 'local.control': ",paste(auxch,collapse=", "),call.=F)
  #  }
  if("adapt" %in% names(local.control)) {
    if(!is.logical(local.control[["adapt"]]) || hasname(local.control[["adapt"]])==F) stop("Component 'adapt' of argument 'local.control' must be a named logical vector",call.=F)
    auxch <- setdiff(names(local.control[["adapt"]]),topG)
    if(length(auxch)>0) stop("Unknown variables in component 'adapt' of argument 'local.control': ",paste(auxch,collapse=", "),call.=F)
    }
  if("max.gestation" %in% names(local.control)) controlCheck("max.gestation",local.control,pset)
  if("min.width" %in% names(local.control)) controlCheck("min.width",local.control,pset)
  if("max.lead" %in% names(local.control)) controlCheck("max.lead",local.control,pset)
  if("sign" %in% names(local.control)) controlCheck("sign",local.control,pset,is.sign=T)
  auxch <- intersect(names(local.control[["min.width"]]),names(local.control[["max.lead"]]))  
  if(length(auxch)>0) {
    for(i in 1:length(auxch)) {
      iminw <- local.control[["min.width"]][[auxch[i]]]
      imaxw <- local.control[["max.lead"]][[auxch[i]]]
      ich <- intersect(names(iminw),names(imaxw))
      if(length(ich)>0) {
        for(j in 1:length(ich)) {
          if(iminw[ich[j]]>imaxw[ich[j]]) stop("Component 'min.width' greater than component 'max.lead' in argument 'local.control': ",auxch[i],"~",ich[j],call.=F)
          }
        }
      }
    }
  nodenam <- c(exogenous,topG)
  if(!is.null(group)) data[,group] <- factor(data[,group])
  xfact <- c()
  for(i in 1:length(nodenam)) {
    if(isQuant(data[,nodenam[i]])==F) {
      xfact <- c(xfact,nodenam[i])
      data[,nodenam[i]] <- factor(data[,nodenam[i]])
      if(length(levels(data[,nodenam[i]]))<2) stop("Variable '",nodenam[i],"' has less than two unique values",sep="",call.=F)
      }
    }
  origdat <- data
  if(log==T) {
    nolog <- xfact
    logtest <- setdiff(nodenam,nolog)            
    if(length(logtest)>0) {
      for(i in 1:length(logtest)) {
        if(sum(data[,logtest[i]]<=0,na.rm=T)>0) {
          nolog <- c(nolog,logtest[i])        
          } else {
          data[,logtest[i]] <- log(data[,logtest[i]])
          }
        }
      }
    if(length(nolog)>0) cat("Logarithm not applied to variables: ",paste(nolog,collapse=", "),sep="","\n")
    }
  difftest <- setdiff(nodenam,xfact)
  if(length(difftest)>0) {
    fine <- ndiff <- 0
    cat("Checking stationarity...")
    flush.console()
    while(fine==0) {
      auxp <- c()
      for(i in 1:length(difftest)) {
       auxdat <- na.omit(data[,difftest[i]])
        ipvl <- unirootTest(difftest[i],group=group,data=data,test=test,combine=combine,k=k,lshort=lshort)$p.value
        ifelse(is.null(ipvl), auxp[i] <- 0, auxp[i] <- ipvl)
        }   
      if(sum(auxp>0.05)>0) {
        if(ndiff<maxdiff) {
          data <- applyDiff(difftest,group=group,data=data,k=rep(1,length(difftest)))                                        
          ndiff <- ndiff+1
          } else {
          fine <- 1
          }
        } else {
        fine <- 1
        }
      }
    if(ndiff==0) {
      cat('\r',"No differentiation performed")
      } else {
      cat('\r',"Order",ndiff,"differentiation performed")    
      }
    cat("\n")
    } else {
    cat('\r',"No differentiation performed")      
    }
  if(is.null(group)) {
    nK <- length(nodenam)+1
    } else {
    nK <- length(nodenam)+length(unique(data[,group]))    
    }
  #for(i in 1:length(nodenam)) {
  #  if(length(na.omit(data[,nodenam[i]]))<nK+1) stop("Imputation impossible: too much missing values for variable '",nodenam[i],"'",call.=F)  #####
  #  }
  auxna <- apply(data[,nodenam],1,function(x){sum(is.na(x))})
  if(sum(auxna)>0) {
    auxOK <- unname(which(auxna<length(nodenam))) 
    if(ndiff>0) {
      if(is.null(group)) {
        auxind <- 1:ndiff
        } else {
        gruppi <- unique(data[,group])
        auxind <- c()
        for(i in 1:length(gruppi)) {
          auxind <- c(auxind,which(data[,group]==gruppi[i])[1:ndiff])          
          }
        }                              
      auxOK <- setdiff(auxOK,auxind)                    
      }
    if(sum(is.na(data[auxOK,]))>0 && maxiter>0) {
      x2imp <- setdiff(nodenam,no.imput)
      if(length(x2imp)>0) {
        data[auxOK,c(group,x2imp)] <- EM.imputation(x=x2imp,group=group,data=data[auxOK,c(group,x2imp),drop=F],tol=tol,maxiter=maxiter)
        #data[auxOK,c(group,nodenam)] <- EM.imputation(x=nodenam,group=group,data=data[auxOK,c(group,nodenam)],tol=tol,maxiter=maxiter)
        }
      }
    }
  nomi <- c()
  optList <- list()
  cat("Start estimation...")
  flush.console()
  for(i in 1:length(model.code)) {
    nomi[i] <- as.character(model.code[[i]])[2]
    auxmess <- paste("Estimating regression model ",i,"/",length(model.code)," (",nomi[i],")",sep="")
    auxdel <- messlen-nchar(auxmess)+1
    if(auxdel>0) {
      cat('\r',auxmess,rep(" ",auxdel))
      } else {
      cat('\r',auxmess)
      }
    flush.console()
    if(is.null(exogenous)) {
      iform <- model.code[[i]]
      } else {
      iform <- formula(paste(paste(as.character(model.code[[i]])[c(2,1,3)],collapse=""),"+",paste(exogenous,collapse="+"),sep=""))
      }                           
    iad <- global.control$adapt
    if("adapt" %in% names(local.control)) {
      if(nomi[i] %in% names(local.control[["adapt"]])) {
        iad <- local.control[["adapt"]][nomi[i]]
        }
      }
    if(length(pset[[nomi[i]]])>0) {
      #if(!is.null(global.control$L)) {
      #  iL <- rep(global.control$L,length(pset[[nomi[i]]]))
      #  names(iL) <- pset[[nomi[i]]]
      #  } else {
      #  iL <- c()
      #  }
      if(!is.null(global.control$max.gestation)) {
        iges <- rep(global.control$max.gestation,length(pset[[nomi[i]]]))
        names(iges) <- pset[[nomi[i]]]
        } else {
        iges <- c()
        }
      if(!is.null(global.control$min.width)) {
        iwd <- rep(global.control$min.width,length(pset[[nomi[i]]]))
        names(iwd) <- pset[[nomi[i]]]
        } else {
        iwd <- c()
        }      
      if(!is.null(global.control$max.lead)) {
        ilead <- rep(global.control$max.lead,length(pset[[nomi[i]]]))       
        names(ilead) <- pset[[nomi[i]]]
        } else {
        ilead <- c()
        }      
      if(!is.null(global.control$sign)) {
        isg <- rep(global.control$sign,length(pset[[nomi[i]]]))
        names(isg) <- pset[[nomi[i]]]
        } else {
        isg <- c()
        }      
      } else {
      iges <- iwd <- ilead <- isg <- c()
      #iL <- c()
      }   
    #if("L" %in% names(local.control)) {
    #  if(nomi[i] %in% names(local.control[["L"]])) {
    #    iauxL <- local.control[["L"]][nomi[i]]
    #    iL[names(iauxL)] <- iauxL
    #    }
    #  }
    if("max.gestation" %in% names(local.control)) {
      if(nomi[i] %in% names(local.control[["max.gestation"]])) {
        iauxges <- local.control[["max.gestation"]][[nomi[i]]]
        iges[names(iauxges)] <- iauxges
        }
      }
    if("min.width" %in% names(local.control)) {
      if(nomi[i] %in% names(local.control[["min.width"]])) {
        iauxwd <- local.control[["min.width"]][[nomi[i]]]
        iwd[names(iauxwd)] <- iauxwd
        }
      }
    if("max.lead" %in% names(local.control)) {
      if(nomi[i] %in% names(local.control[["max.lead"]])) {
        iauxlead <- local.control[["max.lead"]][[nomi[i]]]
        ilead[names(iauxlead)] <- iauxlead
        }
      }
    if("sign" %in% names(local.control)) {
      if(nomi[i] %in% names(local.control[["sign"]])) {
        iauxsg <- local.control[["sign"]][[nomi[i]]]
        isg[names(iauxsg)] <- iauxsg
        }
      }
    optList[[i]] <- list(adapt=iad,max.gestation=iges,min.width=iwd,max.lead=ilead,sign=isg)
    mod0 <- dlaglm(iform,group=group,data=data,adapt=iad,max.gestation=iges,min.width=iwd,max.lead=ilead,sign=isg,selection=global.control$selection)  ### L=iL,
    res[[i]] <- mod0
    messlen <- nchar(auxmess)
    }
  names(optList) <- nomi
  auxmess <- "Estimation completed"
  auxdel <- messlen-nchar(auxmess)+1
  if(auxdel>0) {
    cat('\r',auxmess,rep(" ",auxdel),sep="")
    } else {
    cat('\r',auxmess)
    }
  cat("\n")
  names(res) <- nomi
  finalcode <- list()
  for(i in 1:length(topG)) {
    iscan <- scanForm(res[[topG[i]]]$call$formula)
    icov <- setdiff(iscan$X,c(group,exogenous))
    if(length(icov)==0) icov <- "1"
    finalcode[[i]] <- formula(paste(iscan$y," ~ ",paste(icov,collapse="+"),sep=""))
    res[[topG[i]]]$call$data <- as.name(nameOfData)
    }
  names(finalcode) <- topG
  out <- list(estimate=res[topG],model.code=finalcode,exogenous=exogenous,group=group,log=log,ndiff=ndiff,
    diff.options=diff.options,imput.options=imput.options,selection=global.control$selection,adaptation=optList,
    data.orig=origdat[,c(group,nodenam)],data.used=data[,c(group,nodenam)])
  class(out) <- "dlsem"
  out
  }

# automated plots of lag shapes
auto.lagPlot <- function(x,cumul=FALSE,conf=0.95,plotDir=NULL) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(length(cumul)!=1 || !is.logical(cumul)) stop("Arguent 'cumul' must be a logical value",call.=F)
  if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  if(is.null(plotDir)) plotDir <- getwd()
  for(i in 1:length(x$model.code)) {  
    scan0 <- scanForm(x$model.code[[i]])
    ilagged <- names(scan0$ltype)[which(scan0$ltype!="none")]
    if(length(ilagged)>0) {                                  
      for(j in 1:length(ilagged)) {
        yaux <- rbind(rep(0,3),lagEff(x$estimate[[scan0$y]],ilagged[j],cumul=F,conf=0.95,lag=NULL))
        rownames(yaux)[1] <- "-1"                          
        if(!is.null(yaux)) {
          pdf(file.path(plotDir,paste(scan0$y,"~",ilagged[j],".pdf",sep="")))
          makeShape(yaux,maxlag=NULL,cumul=cumul,conf=conf,ylim=NULL,title=paste(scan0$y," ~ ",ilagged[j],sep=""))
          dev.off()
          }
        }
      }
    }
  cat("Plots saved in ",plotDir,"\n",sep="")
  }

# print method for class dlsem
print.dlsem <- function(x,...) {
  cat("A distributed-lag linear structural equation model","\n")
  n.e <- sum(sapply(edges(makeGraph(x)$graph),length))
  N.e <- sum(sapply(edges(makeGraph(x)$full.graph),length))
  if(!is.null(x$group)) {
    cat(" Group factor: ",x$group,"\n",sep="")
    } else {
    cat(" No group factor","\n")
    }
  if(!is.null(x$exogenous)) {
    cat(" Exogenous variables: ",paste(x$exogenous,collapse=", "),"\n",sep="")
    } else {
    cat(" No exogenous variables","\n")
    }
  if(n.e>0) {
    cat(" ",n.e,"/",N.e," significant edges at 5% level","\n",sep="")
    } else {
    cat(" No significant edges at 5% level","\n")  
    }
  }

# summary method for class dlsem
summary.dlsem <- function(object,...) {
  lapply(object$estimate,summary)
  }

# residuals method for class dlsem
residuals.dlsem <- function(object,...) {
  pred <- list()
  for(i in 1:length(object$estimate)) {
    pred[[i]] <- object$estimate[[i]]$residuals
    }
  names(pred) <- names(object$estimate)
  ind <- lapply(pred,names)
  maxind <- max(as.numeric(unlist(ind)))
  res <- matrix(nrow=maxind,ncol=length(pred))
  colnames(res) <- names(pred)
  rownames(res) <- 1:maxind
  for(i in 1:length(pred)) {
    inam <- names(pred[[i]])
    for(j in 1:length(inam)) {      
      res[as.numeric(inam[j]),i] <- pred[[i]][inam[j]]
      }
    }
  if(!is.null(object$group)) {
    out <- data.frame(object$data.orig[,object$group],res)
    colnames(out)[1] <- object$group
    } else {
    out <- data.frame(res)
    }
  out
  }

# create graph object from dlsem (internal use only)
makeGraph <- function(x,conf=0.95) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  nomi <- names(x$estimate)
  G0 <- G <- new("graphNEL",nodes=nomi,edgemode="directed")
  eSign <- c()
  for(i in 1:length(nomi)) {
    isumm <- summary(x$estimate[[nomi[i]]])$coefficients
    auxnam <- rownames(isumm)
    for(j in 1:length(auxnam)) {
      auxsg <- sign(isumm[auxnam[j],1])
      if(auxnam[j] %in% nomi) {        
        G0 <- addEdge(auxnam[j],nomi[i],G0,1)
        if(isumm[auxnam[j],4]<1-conf) {
          G <- addEdge(auxnam[j],nomi[i],G,1)
          if(auxsg>0) {
            eSign <- c(eSign,"+")
            } else {
            eSign <- c(eSign,"-")
            }
          } else {
          eSign <- c(eSign,"")
          }
        names(eSign)[length(eSign)] <- paste(auxnam[j],"~",nomi[i],sep="")              
        }
      }
    }
  list(graph=G,full.graph=G0,sign=eSign)
  }

# compute the coefficient associated to each edge at different time lags
edgeCoeff <- function(x,lag=NULL,conf=0.95) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  if(!is.null(lag)) {
    for(i in 1:length(lag)) {
      if(!is.numeric(lag[i]) || is.na(lag[i]) || lag[i]<0 || lag[i]!=round(lag[i])) stop("Argument 'lag' must contain non-negative integer numbers only",call.=F)
      }  
    }
  nomi <- names(x$estimate)
  laglen <- c()
  for(i in 1:length(nomi)) {
    isumm <- summary(x$estimate[[nomi[i]]])$coefficients
    auxnam <- rownames(isumm)
    for(j in 1:length(auxnam)) {
      if(auxnam[j] %in% nomi) {
        cumb <- lagEff(x$estimate[[nomi[i]]],x=auxnam[j],cumul=F,conf=conf,lag=NULL)
        laglen <- c(laglen,rownames(cumb)[nrow(cumb)])
        }
      }
    }                                         
  meL <- max(as.numeric(laglen))
  lagOK <- sort(unique(c(lag,0:meL)))
  bList <- vector("list",length=length(lagOK))
  for(i in 1:length(nomi)) {
    isumm <- summary(x$estimate[[nomi[i]]])$coefficients
    auxnam <- rownames(isumm)
    for(j in 1:length(auxnam)) {
      if(auxnam[j] %in% nomi) {
        for(w in 1:length(lagOK)) {
          bList[[w]] <- rbind(bList[[w]],lagEff(x$estimate[[nomi[i]]],x=auxnam[j],cumul=F,conf=conf,lag=lagOK[w]))
          rownames(bList[[w]])[nrow(bList[[w]])] <- paste(auxnam[j],"~",nomi[i],sep="")
          }
        }
      }
    if(!is.null(bList)) names(bList) <- lagOK
    }
  for(i in 1:length(bList)) {
    auxnam <- rownames(bList[[i]])
    newnam <- c()
    for(j in 1:length(auxnam)) {
      newnam[j] <- paste(rev(strsplit(auxnam[j],"~")[[1]]),collapse="~")
      }
    rownames(bList[[i]]) <- newnam
    }
  for(i in 1:length(bList)) {
    colnames(bList[[i]]) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
    }
  if(is.null(lag)) {
    bList
    } else {
    voidB <- bList[[1]]
    voidB[] <- 0
    bList2 <- vector("list",length=length(lag))
    names(bList2) <- lag
    for(i in 1:length(bList2)) {
      if(lag[i]<=max(lagOK)) {
        bList2[[i]] <- bList[[as.character(lag[i])]]
        } else {
        bList2[[i]] <- voidB
        }
      }
    bList2
    }
  }

# plot method for class dlsem
plot.dlsem <- function(x,conf=0.95,style=2,node.col=NULL,font.col=NULL,border.col=NULL,edge.col=NULL,edge.lab=NULL,...) {
  if((style %in% c(0,1,2))==F) stop("Argument 'style' must be either '0' (plain), '1' (significance shown), or '2' (signs shown)",call.=F)
  defpar <- par()[setdiff(names(par()),c("cin","cra","csi","cxy","din","page"))]
  G <- makeGraph(x,conf=conf)
  nomi <- nodes(G$full.graph)  
  #####
  cutString <- function(x,l) {
    n <- nchar(x)
    k <- ceiling(n/l)
    res <- c()
    for(i in 1:k) {
      res[i] <- substr(x,1+(i-1)*l,i*l)
      }
    paste(res,collapse="\n")
    }
  #####
  nAttr <- list()
  nAttr$shape <- rep("ellipse",length(nomi))
  nAttr$fontsize <- rep(14,length(nomi))
  nAttr$height <- rep(2.5,length(nomi))
  nAttr$width <- rep(4,length(nomi))
  nAttr$label <- sapply(nomi,cutString,l=12)  ##### maximum number of characters: 12
  for(i in 1:length(nAttr)) {
    names(nAttr[[i]]) <- nomi
    }
  if(!is.null(node.col)) nAttr$fillcolor <- node.col
  if(!is.null(font.col)) nAttr$fontcolor <- font.col
  if(!is.null(border.col)) nAttr$color <- border.col
  eAttr <- list()                                                     
  if(!is.null(edge.lab)) {
    eAttr$label <- edge.lab
    #eAttr$labelfontsize <- rep(14,length(edge.lab))
    #names(eAttr$labelfontsize) <- names(edge.lab)
    }
  if(!is.null(edge.col)) {
    eAttr$color <- edge.col     
    eAttr$color[which(G$sign=="")] <- NA
    } else {
    eCol <- G$sign
    if(style==1) {
      eCol[which(G$sign=="+")] <- "grey30"
      eCol[which(G$sign=="-")] <- "grey30"
      eCol[which(G$sign=="")] <- "grey75"
      } else if(style==2) {
      eCol[which(G$sign=="+")] <- "green4"
      eCol[which(G$sign=="-")] <- "tomato3"
      eCol[which(G$sign=="")] <- "grey75"
      } else {
      eCol[] <- "grey30"
      }
    eAttr[[length(eAttr)+1]] <- eCol
    }
  names(eAttr)[length(eAttr)] <- "color"
  #eStyl <- G$sign
  #eStyl[which(G$sign=="+")] <- "solid"
  #eStyl[which(G$sign=="-")] <- "dashed"
  #eAttr[[length(eAttr)+1]] <- eStyl
  #names(eAttr)[length(eAttr)] <- "style"  
  plot(G$full.graph,"dot",nodeAttrs=nAttr,edgeAttrs=eAttr,attrs=list(edge=list(color="grey25",arrowsize=0.4)),...)
  par(defpar)
  }   
      
# parent sets (internal use only)
edgeIn <- function(G) {
  xedg <- edges(G)
  pset <- vector("list",length=length(xedg))
  nomi <- names(pset) <- names(xedg)
  for(i in 1:length(nomi)) pset[[nomi[i]]] <- character(0)
  for(i in 1:length(nomi)) {
    for(j in 1:length(nomi)) {
      if(nomi[i] %in% xedg[[nomi[j]]]) pset[[nomi[i]]] <- c(pset[[nomi[i]]],nomi[j])
      }
    }
  pset
  }

# node parents (internal use only)
nodeParen <- function(nodeName,G) {
  xedg <- edges(G)
  pset <- character(0)
  nomi <- names(xedg)
  for(j in 1:length(nomi)) {
    if(nodeName %in% xedg[[nomi[j]]]) pset <- c(pset,nomi[j])
    }
  pset
  }

# node children (internal use only)
nodeChld <- function(nodeName,G) {
  edges(G)[[nodeName]]
  }

# node markov blanket (internal use only)
nodeMB <- function(nodeName,G) {
  auxch <- nodeChld(nodeName,G)
  auxpar <- nodeParen(nodeName,G)
  auxp <- c()
  if(length(auxch)>0) {
    for(i in 1:length(auxch)) {
      auxp <- c(auxp,nodeParen(auxch[i],G))
      }
    }
  setdiff(unique(c(auxpar,auxch,auxp)),nodeName)
  }

# node ancestors (internal use only)
nodeAnces <- function(nodeName,G) {
  eList <- edgeIn(G)
  xpar <- aux.xpar <- eList[[nodeName]]
  while(length(aux.xpar)>0) {
    newpar <- c()
    for(i in 1:length(aux.xpar)) {
      newpar <- c(newpar,eList[[aux.xpar[i]]])
      }
    xpar <- unique(c(xpar,newpar))
    aux.xpar <- newpar
    }
  unique(xpar)
  }

# node descendants (internal use only)
nodeDescen <- function(nodeName,G) {
  eList <- edges(G)
  xpar <- aux.xpar <- eList[[nodeName]]
  while(length(aux.xpar)>0) {
    newpar <- c()
    for(i in 1:length(aux.xpar)) {
      newpar <- c(newpar,eList[[aux.xpar[i]]])
      }
    xpar <- c(xpar,newpar)
    aux.xpar <- newpar
    }
  unique(xpar)
  }

# find topological order (internal use only)
topOrder <- function(G) {
  parSet <- edgeIn(G)
  nomi <- names(parSet)
  L <- c()
  S <- nomi[which(sapply(parSet,length)==0)] 
  while(length(S)>0) {
    xaux <- S[1]
    S <- setdiff(S,xaux)
    L <- c(L,xaux)
    sch <- c()
    for(j in 1:length(parSet)) {
      if(xaux %in% parSet[[j]]) sch <- c(sch,nomi[j])
      }
    if(length(sch)>0) {
      for(j in 1:length(sch)) {
        parSet[[sch[j]]] <- setdiff(parSet[[sch[j]]],xaux)
        if(length(parSet[[sch[j]]])==0) S <- c(S,sch[j])  
        }
      }
    }
  if(sum(sapply(parSet,length))==0) L else NULL
  }

# check conditional independence
isIndep <- function(x,var1=NULL,var2=NULL,given=NULL,conf=0.95,use.ns=FALSE) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(is.null(var1) || is.na(var1)) stop("Argument 'var1' is missing",call.=F)
  if(!is.null(var1) && length(var1)!=1) stop("Argument 'var1' must have length 1",call.=F)
  if(is.null(var2) || is.na(var2)) stop("Argument 'var2' is missing",call.=F)  
  if(!is.null(var2) && length(var2)!=1) stop("Argument 'var2' must have length 1",call.=F)
  if(!is.null(given) && is.na(given)) stop("Argument 'given' is missing",call.=F)  
  if(length(use.ns)!=1 || !is.logical(use.ns)) stop("Argument 'use.ns' must be a logical value",call.=F)
  Gobj <- makeGraph(x,conf=conf)
  if(use.ns==F) {
    G <- Gobj$graph
    } else {
    G <- Gobj$full.graph
    }
  nomi <- nodes(G)
  auxcheck <- setdiff(c(var1,var2,given),nomi)
  if(length(auxcheck)>0) stop("Unknown variable '",auxcheck[1],"'",sep="",call.=F)
  if(length(unlist(edges(G)))>0) {
    Gm <- moralize(ancestralGraph(c(var1,var2,given),G))
    pset <- edges(Gm)                      
    xedg <- c()
    for(i in 1:length(pset)) {
      if(length(pset[[i]])>0) {
        for(j in 1:length(pset[[i]])) {
          xedg <- rbind(xedg,c(names(pset)[i],pset[[i]][j])) 
          }
        }
      }                              
    if(!is.null(xedg)) {
      xedg <- rbind(xedg,xedg[,2:1])
      nomi <- sort(unique(xedg))
      borde <- var1
      reached <- c()
      neighb <- function(x,node) {
        Ne <- c(x[which(x[,1]==node),2],x[which(x[,2]==node),1])
        sort(unique(Ne))
        }                     
      while(length(borde)>0) {
        reached <- c(reached,borde)
        fan_borde <- c()
        for(i in 1:length(borde)) {                      
          auxne <- neighb(xedg,borde[i])                 
          fan_borde <- c(fan_borde,auxne) 
          }                              
        borde <- setdiff(fan_borde,c(reached,given))   
        if(length(intersect(borde,var2))>0) break()
        }
      ifelse(length(borde)>0, res <- F, res <- T)
      } else {
      res <- T
      }
    } else {
    res <- T
    }
  res
  }

# matrix of edges (internal use only)
edgeMat <- function(x,conf,full) {
  if(full==T) {
    xedg <- edges(makeGraph(x,conf=conf)$full.graph)
    } else {
    xedg <- edges(makeGraph(x,conf=conf)$graph)
    }
  res <- c()
  if(length(xedg)>0) {
    for(i in 1:length(xedg)) {
      if(length(xedg[[i]])>0) {
        res <- rbind(res,cbind(rep(names(xedg)[i],length(xedg[[i]])),xedg[[i]]))
        }
      }
    }
  res
  }

# find directed path (internal use only)
dpathFind <- function(G,from,to) {
  auxedg <- edges(G)
  if(to %in% nodeDescen(from,G)) {
    auxmat <- c()
    for(i in 1:length(auxedg)) {
      if(length(auxedg[[i]])>0) {
        for(j in 1:length(auxedg[[i]])) {
          auxmat <- rbind(auxmat,c(names(auxedg)[i],auxedg[[i]][j]))
          }
        }
      }
    auxchld <- auxmat[which(auxmat[,1]==from),2]
    pathList <- list()
    for(i in 1:length(auxchld)) pathList[[i]] <- c(from,auxchld[i])
    endcheck <- function(x) {
      res <- rep(0,length(x))                 
      for(i in 1:length(x)) {                
        if(rev(x[[i]])[1]==to) res[i] <- 1
        }
      res
      }                        
    isOK <- endcheck(pathList)               
    while(sum(isOK)<length(isOK)) {
      auxind <- which(isOK==0)[1]         
      auxchld <- auxmat[which(auxmat[,1]==rev(pathList[[auxind]])[1]),2]
      auxpath <- pathList[[auxind]]
      if(length(auxchld)>0) {
        pathList[[auxind]] <- c(auxpath,auxchld[1])
        if(length(auxchld)>1) {
          for(i in 2:length(auxchld)) {
            pathList[[length(pathList)+1]] <- c(auxpath,auxchld[i])
            }
          }   
        } else {
        pathList <- pathList[-auxind]
        }                       
      isOK <- endcheck(pathList)  
      }
    pathList
    } else {
    NULL
    }
  } 

# find lag to sum (internal use only)
findlag2sum <- function(x,lag) {
  g <- length(x)                             
  out <- list()
  for(w in 1:length(lag)) {
    mycode <- "res <- c(); "
    for(i in 1:g) {
      mycode <- paste(mycode,"for(k",i," in 1:length(x[[",i,"]])) { ; ",sep="")
      }
    mycode <- paste(mycode,paste("xaux <- c(",paste("x[[",1:g,"]][k",1:g,"]",collapse=",",sep=""),"); ",sep=""),sep="")
    mycode <- paste(mycode,"if(sum(xaux)==lag[w]) { res <- rbind(res,xaux) }; ",sep="")
    mycode <- paste(mycode,paste(rep("}",length(x)),collapse="; "),sep="")
    eval(parse(text=mycode))
    if(is.null(res)) {
      res <- matrix(nrow=0,ncol=g)
      } else {
      rownames(res) <- NULL
      }
    colnames(res) <- names(x)             
    out[[w]] <- res
    }
  names(out) <- lag
  out
  }

# computation of causal effects
causalEff <- function(x,from=NULL,to=NULL,lag=NULL,cumul=FALSE,conf=0.95,use.ns=FALSE) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(is.null(from) || is.na(from)) stop("Argument 'from' is missing",call.=F)
  if(is.null(to) || is.na(to)) stop("Argument 'to' is missing",call.=F)
  if(!is.character(from)) stop("Invalid argument 'from'",call.=F)
  if(!is.character(to) || length(to)!=1) stop("Invalid argument 'to'",call.=F)
  if(length(cumul)!=1 || !is.logical(cumul)) stop("Arguent 'cumul' must be a logical value",call.=F)
  if(length(use.ns)!=1 || !is.logical(use.ns)) stop("Argument 'use.ns' must be a logical value",call.=F)
  if(!is.null(lag)) {
    for(i in 1:length(lag)) {
      if(!is.numeric(lag[i]) || is.na(lag[i]) || lag[i]<0 || lag[i]!=round(lag[i])) stop("Argument 'lag' must contain non-negative integer numbers only",call.=F)
      }
    }
  auxcheck <- setdiff(c(from,to),names(x$estimate))
  if(length(auxcheck)>0) {
    auxcntx <- intersect(auxcheck,x$exogenous)
    if(length(auxcntx)>0) {
      stop("Variable '",auxcntx[1],"' is exogenous and cannot appear in argument 'from' or 'to'",call.=F)
      } else {
      stop("Unknown variable '",auxcheck[1],"'",sep="",call.=F)
      }
    }
  Gobj <- makeGraph(x,conf=conf)
  if(use.ns==F) {
    G <- Gobj$graph
    } else {
    G <- Gobj$full.graph
    }
  nomi <- nodes(G)
  if(length(setdiff(from,nomi))>0) stop("Unknown variable '",setdiff(from,nomi)[1],"'",sep="",call.=F)
  if((to %in% nomi)==F) stop("Unknown variable '",to,"'",sep="",call.=F)
  isOK <- rep(1,length(from))
  for(i in 1:length(from)) {
    if((to %in% nodeDescen(from[i],G))==F) isOK[i] <- 0
    }
  if(sum(isOK)>0) {
    from <- from[which(isOK==1)]
    mycol1 <- mycol2 <- rep(NA,length(nomi))
    names(mycol1) <- names(mycol2) <- nomi
    nodemed <- c()
    pathList <- vector("list",length=length(from))
    for(i in 1:length(from)) {
      pathList[[i]] <- dpathFind(G,from=from[i],to=to)
      nodemed <- c(nodemed,setdiff(unlist(pathList[[i]]),c(from,to)))
      }
    nodemed <- unique(nodemed)
    names(pathList) <- from
    auxdel <- which(sapply(pathList,is.null)==T)
    if(length(auxdel)>0) {
      pathList <- pathList[-auxdel]
      from <- from[-auxdel]
      }
    pset <- edgeIn(G)
    nodecond <- setdiff(unlist(pset[c(to,nodemed)]),c(from,nodemed))
    nodebarr <- setdiff(nomi,c(from,to,nodemed,nodecond))
    mycol1[c(from,to)] <- mycol2[c(from,to)] <- "grey20"
    mycol1[nodemed] <- mycol2[nodemed] <- "grey20"
    mycol1[nodecond] <- "navy"
    mycol2[nodecond] <- "grey70"
    mycol1[nodebarr] <- mycol2[nodebarr] <- "grey70"
    xedg <- edges(G)
    ednam <- list()
    for(i in 1:length(xedg)) {
      if(length(xedg[[i]])>0) ednam[[i]] <- paste(names(xedg)[i],"~",xedg[[i]],sep="")
      }
    ednam <- unlist(ednam)
    eddel <- c()
    for(i in 1:length(nodebarr)) {
      eddel <- c(eddel,paste(nodebarr[i],"~",setdiff(nomi,nodebarr[i]),sep=""),paste(setdiff(nomi,nodebarr[i]),"~",nodebarr[i],sep=""))
      }
    for(i in 1:length(nodecond)) {
      eddel <- c(eddel,paste(nodecond[i],"~",setdiff(nomi,nodecond[i]),sep=""),paste(setdiff(nomi,nodecond[i]),"~",nodecond[i],sep=""))
      }
    edcol <- rep("grey70",length(Gobj$sign))
    names(edcol) <- names(Gobj$sign)
    edcol[intersect(setdiff(ednam,eddel),names(which(Gobj$sign=="+")))] <- "green4"
    edcol[intersect(setdiff(ednam,eddel),names(which(Gobj$sign=="-")))] <- "tomato3"
    newPathList <- list()
    for(i in 1:length(pathList)) {
      newPathList <- c(newPathList,pathList[[i]])
      }                                                        
    laglen <- list()
    for(i in 1:length(newPathList)) {
      jlaglen <- list()
      for(j in 2:length(newPathList[[i]])) {
        auxnam <- paste(newPathList[[i]][j],"~",newPathList[[i]][j-1],sep="")
        auxeff <- lagEff(x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf,lag=NULL)
        auxpos <- which(auxeff[,1]!=0)
        if(length(auxpos)>0) {
          jlaglen[[j-1]] <- as.numeric(rownames(auxeff)[auxpos])
          } else {
          jlaglen[[j-1]] <- NA
          }
        }
      laglen[[i]] <- c(min(jlaglen[[1]],na.rm=T),sum(sapply(jlaglen,max,na.rm=T)))
      }
    meL <- max(unlist(laglen))+1        
    lagOK <- 0:meL
    mycoeff <- edgeCoeff(x,lag=lagOK,conf=conf)         
    sd_calc <- function(muvet,sdvet) { sqrt(prod(muvet^2+sdvet^2)-prod(muvet^2)) }
    quan <- -qnorm((1-conf)/2)
    bhat <- list()
    for(i in 1:length(mycoeff)) {
      bhat[[i]] <- matrix(nrow=nrow(mycoeff[[i]]),ncol=2)
      for(j in 1:nrow(mycoeff[[i]])) {
        auxeval <- mycoeff[[i]][j,1]
        auxsd <- (mycoeff[[i]][j,3]-mycoeff[[i]][j,1])/quan
        bhat[[i]][j,] <- c(auxeval,auxsd)
        }
      rownames(bhat[[i]]) <- rownames(mycoeff[[i]])
      }
    names(bhat) <- names(mycoeff)
    outList <- list()                       
    for(i in 1:length(newPathList)) {           
      outList[[i]] <- matrix(nrow=length(lagOK),ncol=3)
      rownames(outList[[i]]) <- lagOK
      colnames(outList[[i]]) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
      auxbetalag <- list()
      for(j in 2:length(newPathList[[i]])) {
        auxnam <- paste(newPathList[[i]][j],"~",newPathList[[i]][j-1],sep="")
        auxeff <- lagEff(x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf,lag=lagOK)
        auxpos <- which(auxeff[,1]!=0)
        if(length(auxpos)>0) {                                 
          auxbetalag[[j-1]] <- as.numeric(rownames(auxeff)[auxpos])
          } else {
          auxbetalag[[j-1]] <- 0
          }
        names(auxbetalag)[j-1] <- auxnam
        }
      lagsumMat <- findlag2sum(auxbetalag,lagOK)
      for(j in 1:length(lagsumMat)) {
        auxlag <- as.character(lagOK[j])                           
        auxind <- lagsumMat[[j]]                          
        if(nrow(auxind)>0) {
          auxres <- array(dim=c(nrow(auxind),ncol(auxind),2))
          for(w1 in 1:nrow(auxind)) {                          
            for(w2 in 1:ncol(auxind)) {
              auxres[w1,w2,1] <- bhat[[as.character(auxind[w1,w2])]][colnames(auxind)[w2],1]
              auxres[w1,w2,2] <- bhat[[as.character(auxind[w1,w2])]][colnames(auxind)[w2],2]
              }
            }                                                           
          muprod <- sdprod <- c()
          for(w in 1:nrow(auxind)) {                               
            muprod[w] <- prod(auxres[w,,1])                              
            sdprod[w] <- sd_calc(auxres[w,,1],auxres[w,,2])
            }
          mupath <- sum(muprod)
          sdpath <- sqrt(sum(sdprod^2))
          outList[[i]][j,] <- c(mupath,mupath-quan*sdpath,mupath+quan*sdpath)
          } else {
          outList[[i]][j,] <- rep(0,3)        
          }
        }
      }
    names(outList) <- sapply(newPathList,function(x){paste(x,collapse="*")})
    out <- matrix(nrow=length(lagOK),ncol=3)
    rownames(out) <- lagOK
    colnames(out) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
    for(j in 1:length(lagOK)) {
      auxover <- rep(0,3)
      for(i in 1:length(outList)) {
        auxover <- auxover+outList[[i]][j,]
        }
      out[j,] <- auxover
      }
    outList[[length(outList)+1]] <- out
    names(outList)[[length(outList)]] <- "overall"
    if(cumul==T) {
      for(i in 1:length(outList)) {
        if(nrow(outList[[i]])>1) {
          for(j in 2:nrow(outList[[i]])) {        
            outList[[i]][j,] <- outList[[i]][j-1,]+outList[[i]][j,]
            }
          }
        }    
      }
    if(is.null(lag)) {
      outList
      } else {
      lagSel <- lag
      auxOlag <- which(lag>max(lagOK))
      if(length(auxOlag)>0) lagSel[auxOlag] <- max(lagOK)  
      outList2 <- outList
      for(i in 1:length(outList)) {
        outList2[[i]] <- outList[[i]][as.character(lagSel),,drop=F]
        rownames(outList2[[i]]) <- lag
        }
      outList2
      }
    } else {
    auxanc <- nodeAnces(to,makeGraph(x)$full.graph)
    if(length(intersect(from,auxanc))>0) {
      stop("No paths found connecting the selected variables. Try to reduce 'conf' or to set 'use.ns' to TRUE",call.=F)
      } else {
      stop("No paths exist connecting the selected variables",call.=F)
      }
    }
  }

# fit indices
modelFit <- function(x) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  xfit <- x$estimate     
  deg <- n <- expdev <- resdev <- Rsq <- ll <- aic <- bic <- mdl <- c()
  for(i in 1:length(xfit)) {
    resdev[i] <- deviance(xfit[[i]])
    auxss <- anova(xfit[[i]])$'Sum Sq'
    expdev[i] <- sum(auxss[1:(length(auxss)-1)])
    n[i] <- length(residuals(xfit[[i]]))
    deg[i] <- xfit[[i]]$df.residual
    Rsq[i] <- summary(xfit[[i]])$'r.squared'
    ll[i] <- modscore(xfit[[i]],"ll")
    aic[i] <- modscore(xfit[[i]],"aic")
    bic[i] <- modscore(xfit[[i]],"bic")
    mdl[i] <- modscore(xfit[[i]],"mdl")
    }  
  out <- c(Rsq,sum(expdev)/(sum(expdev)+sum(resdev)))
  names(out) <- c(names(xfit),"(overall)")
  out0 <- c(ll,sum(ll))
  names(out0) <- c(names(xfit),"(overall)")
  out1 <- c(aic,sum(aic))
  names(out1) <- c(names(xfit),"(overall)")
  out2 <- c(bic,sum(bic))
  names(out2) <- c(names(xfit),"(overall)")
  out3 <- c(mdl,sum(mdl))
  names(out3) <- c(names(xfit),"(overall)")
  list('Rsq'=out,'logL'=out0,'AIC'=out1,'BIC'=out2,'MDL'=out3)
  }

