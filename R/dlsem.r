# lag function (internal use only)
lagFun <- function(x,maxlag,addNA=T) {
  out <- x
  if(maxlag>0) {
    for(w in 1:maxlag) {
      if(addNA==T) {
        z <- rep(NA,length(x))
        } else {
        z <- rep(0,length(x))        
        }
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
  }

# Z matrix (internal use only)
Zmat <- function(x,type,theta) {
  if(type=="quec") {
    minlag <- theta[1]
    maxlag <- theta[2]
    H <- c()
    for(i in 0:maxlag) {
      if(i<minlag | i>maxlag) {
        H[i+1] <- 0
        } else {
        H[i+1] <- -4/(maxlag-minlag+2)^2*(i^2-(minlag+maxlag)*i+(minlag-1)*(maxlag+1))
        }
      }
    lagFun(x,maxlag)%*%matrix(H,ncol=1)
    } else if(type=="qdec") {
    minlag <- theta[1]
    maxlag <- theta[2]
    H <- c()
    for(i in 0:maxlag) {
      if(i<minlag | i>maxlag) {
        H[i+1] <- 0
        } else {
        H[i+1] <- (i^2-2*maxlag*i+maxlag^2)/(maxlag-minlag)^2
        }
      }
    lagFun(x,maxlag)%*%matrix(H,ncol=1)
    #} else if(type=="trap") {
    #xa <- theta[1]
    #xb <- theta[2]
    #xc <- theta[3]
    #xd <- theta[4]
    #H <- c()
    #for(i in 0:xd) {
    #  if(i>=xa & i<=xb) {
    #    H[i+1] <- (i-xa+1)/(xb-xa+1)
    #    } else if(i>xb & i<=xc) {
    #    H[i+1] <- 1
    #    } else if(i>=xc & i<=xd) {
    #    H[i+1] <- (i-xd-1)/(xc-xd-1)
    #    } else {
    #    H[i+1] <- 0
    #    }
    #  }
    #lagFun(x,xd)%*%matrix(H,ncol=1)    
    } else if(type=="gamma") {
    delta <- theta[1]
    lambda <- theta[2]
    xa <- 0
    xb <- qgamma(0.99,1/(1-delta),-log(lambda))+xa
    H <- c()
    for(i in 0:xb) {
      if(i>=xa) {
        bnum <- (i-xa+1)^(delta/(1-delta))*lambda^(i-xa)
        xM <- (delta/(delta-1))/log(lambda)+xa-1
        bden <- (xM-xa+1)^(delta/(1-delta))*lambda^(xM-xa)
        H[i+1] <- bnum/bden
        } else {
        H[i+1] <- 0     
        }
      }                             
    lagFun(x,xb)%*%matrix(H,ncol=1)
    }
  }

# compute aic (internal use only)
modFit <- function(x,method="aic") {
  res <- residuals(x)
  yfit <- fitted(x)
  n <- length(res)
  s2 <- summary(x)$sigma^2
  ll <- -0.5*n*log(2*pi*s2)-0.5/s2*sum(res^2)
  npar <- 2*(length(x$coefficients)+1)
  if(method=="aic") {
    -2*ll+2*npar
    } else if(method=="bic") {
    -2*ll+log(n)
    } else if(method=="mdl") {
    (n-npar)*log(sum(res^2)/npar)+npar*log(sum(yfit^2))+(n-npar-1)*log(n/(n-npar))-(npar+1)*log(npar)
    } else {
    stop("Invalid method. Chose one among 'aic', 'bic' and 'mdl'",call.=F) ####
    }
  }

# least squares estimate (internal use only)
dlagLS <- function(y,X,group,data,type,theta,L) {
  panelList <- list()
  if(is.null(group)) {
    xtime <- panelList[[1]] <- 1:nrow(data)
    Z <- data.frame(rep(NA,nrow(data)))
    auxcall <- c()
    Znam <- y
    if(length(X)>0) {
      for(i in 1:length(X)) {
        if(type[i]=="quec") {
          iZ <- Zmat(data[,X[i]],type[i],theta[[i]])
          Znam <- c(Znam,paste("theta0_quec.",X[i],sep=""))
          } else if(type[i]=="qdec") {
          iZ <- Zmat(data[,X[i]],type[i],theta[[i]])
          Znam <- c(Znam,paste("theta0_qdec.",X[i],sep=""))          
          #} else if(type[i]=="trap") {
          #iZ <- Zmat(data[,X[i]],type[i],theta[[i]])
          #Znam <- c(Znam,paste("theta0_trap.",X[i],sep=""))
          } else if(type[i]=="gamma") {
          iZ <- Zmat(data[,X[i]],type[i],theta[[i]])
          Znam <- c(Znam,paste("theta0_gamma.",X[i],sep=""))
          } else {
          iZ <- data[,X[i]]
          Znam <- c(Znam,X[i])
          }
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
    if(nrow(na.omit(Z))==0) stop("Insufficient number of observations",call.=F)
    mod0 <- lm(formula(form0),data=Z)
    } else {
    allZ <- xtime <- c()
    gruppi <- unique(data[,group])
    for(w in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[w])
      xtime[auxind] <- 1:length(auxind)
      panelList[[w]] <- auxind
      Z <- data.frame(rep(NA,length(auxind)))
      auxcall <- c()
      Znam <- y
      if(length(X)>0) {
        for(i in 1:length(X)) {
          if(type[i]=="quec") {
            iZ <- Zmat(data[auxind,X[i]],type[i],theta[[i]])
            Znam <- c(Znam,paste("theta0_quec.",X[i],sep=""))
            } else if(type[i]=="qdec") {
            iZ <- Zmat(data[auxind,X[i]],type[i],theta[[i]])
            Znam <- c(Znam,paste("theta0_qdec.",X[i],sep=""))
            #} else if(type[i]=="trap") {
            #iZ <- Zmat(data[auxind,X[i]],type[i],theta[[i]])
            #Znam <- c(Znam,paste("theta0_trap.",X[i],sep=""))
            } else if(type[i]=="gamma") {
            iZ <- Zmat(data[auxind,X[i]],type[i],theta[[i]])
            Znam <- c(Znam,paste("theta0_gamma.",X[i],sep=""))
            } else {
            iZ <- data[auxind,X[i]]
            Znam <- c(Znam,X[i])
            }
          Z <- cbind(Z,iZ)
          }
        }
      colnames(Z) <- Znam
      Z[,1] <- data[auxind,y]
      allZ <- rbind(allZ,Z)
      }
    allZ <- cbind(data[,group],allZ)
    colnames(allZ)[1] <- group
    if(nrow(na.omit(allZ))==0) stop("Insufficient number of observations",call.=F)                                        
    names(panelList) <- gruppi
    if(length(X)>0) {
      form0 <- paste(y," ~ -1+factor(",group,")+",paste(colnames(allZ)[-c(1:2)],collapse="+"),sep="")
      } else {
      form0 <- paste(y," ~ -1+factor(",group,")",sep="")
      }
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
  if(L>0) {
    res0 <- residuals(mod0)
    W <- matrix(nrow=length(res0),ncol=length(res0))
    for(i in 1:length(res0)) {
      for(j in 1:length(res0)) {
        ijlag <- abs(xtime[as.numeric(names(res0)[i])]-xtime[as.numeric(names(res0)[j])])
        W[i,j] <- max(0,(1-ijlag/(1+L)))*res0[i]*res0[j]
        }
      }
    Xmat <- model.matrix(mod0)
    Imat <- solve(t(Xmat)%*%Xmat)
    mod0$vcov <- Imat%*%t(Xmat)%*%W%*%Xmat%*%Imat
    } else {
    mod0$vcov <- bvar*summary(mod0)$sigma^2
    }
  mod0$panel <- panelList
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
        if(grepl("quec\\(",auX[i])) {                           
          istr <- gsub("\\)","",gsub("quec\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "quec"
          lpar[[i]] <- as.numeric(istr[-1])
          } else if(grepl("qdec\\(",auX[i])) {                           
          istr <- gsub("\\)","",gsub("qdec\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "qdec"
          lpar[[i]] <- as.numeric(istr[-1])   
          #} else if(grepl("trap\\(",auX[i])) {
          #istr <- gsub("\\)","",gsub("trap\\(","",strsplit(auX[i],",")[[1]]))
          #lnames[i] <- istr[1]
          #ltype[i] <- "trap"
          #lpar[[i]] <- as.numeric(istr[-1])            
          } else if(grepl("gamma\\(",auX[i])) {
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

# gamma constrains (internal use only)
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
searchGrid <- function(maxgs,minwd,maxwd,laglim,lag.type) {
  if(lag.type=="gamma") {
    auxind <- which(gammaCons$sx<=maxgs & gammaCons$dx<=laglim & gammaCons$dx-gammaCons$sx>=minwd & gammaCons$dx-gammaCons$sx<=maxwd,arr.ind=T)
    n <- nrow(auxind)
    if(n>0) {
      res <- c()
      for(i in 1:nrow(auxind)) {
        res <- rbind(res,as.numeric(c(rownames(gammaCons$sx)[auxind[i,1]],colnames(gammaCons$sx)[auxind[i,2]])))
        }
      outind <- seq(1,nrow(auxind),by=ceiling(n/500))
      res[outind,]
      }
    #} else if(lag.type=="trap") {
    #####
    } else {
    auxmat <- c()
    for(i in 0:laglim) {
      for(j in i:laglim) { 
        auxmat <- rbind(auxmat,c(i,j))
        }   
      }
    auxind <- which(auxmat[,1]<=maxgs & auxmat[,2]<=laglim & auxmat[,2]-auxmat[,1]>=minwd & auxmat[,2]-auxmat[,1]<=maxwd)
    n <- length(auxind)
    if(n>0) {
      auxmat[auxind,]
      }
    }
  }
  
# fit a distributed-lag regression model (internal use only)
dlaglm <- function(formula,group=NULL,data,log=FALSE,L=0,adapt=FALSE,max.gestation=NULL,min.width=NULL,max.width=NULL,sign=NULL,selection="aic") {
  auxscan <- scanForm(formula)
  y <- auxscan$y
  lagPar <- auxscan$lpar
  lagType <- auxscan$ltype
  lagNam <- names(lagPar)  
  if(length(lagNam)>0) {
    if(!is.null(group) && (group %in% lagNam)) stop("Variable ",group," is defined as both a group factor and a covariate",call.=F) 
    auxcheck <- setdiff(lagNam,colnames(data))
    if(length(auxcheck)>0) stop("Unknown variable: ",auxcheck[1],call.=F)    
    auxdupl <- duplicated(lagNam)
    if(sum(auxdupl)>0) {
      stop("Duplicated covariates: ",paste(lagNam[auxdupl],collapse=", "),call.=F)
      }
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
  if(log==T) {
    nodenam <- c(y,lagNam)
    nolog <- c()
    for(i in 1:length(nodenam)) {
      if(sum(data[,nodenam[i]]<=0,na.rm=T)==0) {
        data[,nodenam[i]] <- log(data[,nodenam[i]])
        } else {
        nolog <- c(nolog,nodenam[i])
        }
      }                                                                                             
    if(length(nolog)>0) warning("Log not applied to variables: ",paste(nolog,collapse=", "),call.=F)
    }
  if(adapt==F) {
    modOK <- dlagLS(y=y,X=lagNam,group=group,data=data,type=lagType,theta=lagPar,L=L)
    if(length(c(group,auxscan$X))>0) {
      modOK$call <- paste(y," ~ ",paste(c(group,auxscan$X),collapse="+"),sep="")
      } else {
      modOK$call <- paste(y," ~ 1",sep="")
      }
    } else {
    if(is.null(group)) {
      laglimit <- max(0,nrow(na.omit(data))-3)
      } else {
      auxlaglim <- c()
      for(i in 1:length(gruppi)) {
        auxlaglim <- nrow(na.omit(data[which(data[,group]==gruppi[i]),]))
        }
      laglimit <- max(0,min(auxlaglim,na.rm=T)-3)
      }
    bestPar <- vector("list",length=length(lagNam))
    names(bestPar) <- lagNam
    xOK <- c()
    fine <- 0
    while(fine==0) {
      xtest <- setdiff(lagNam,xOK)                 
      if(length(xtest)>0) {
        currentAIC <- currentP <- rep(NA,length(xtest))
        names(currentAIC) <- names(currentP) <- xtest
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
            if(xtest[i] %in% names(max.width)) {
              maxwd <- max.width[xtest[i]]
              } else {
              maxwd <- laglimit
              }
            auxcons <- searchGrid(maxgs,minwd,maxwd,laglimit,lagType[xtest[i]])
            aic0 <- bhat0 <- pval0 <- c()               
            for(j in 1:nrow(auxcons)) {
              testType <- lagType[c(xOK,xtest[i])]     
              testPar <- lagPar[c(xOK,xtest[i])]                       
              testPar[[xtest[i]]] <- auxcons[j,]      
              mod0 <- dlagLS(y=y,X=names(testPar),group=group,data=data,type=testType,theta=testPar,L=0)
              summ0 <- summary(mod0)$coefficients
              ixstr <- paste("theta0_",lagType[xtest[i]],".",xtest[i],sep="")                                                        
              if(ixstr %in% rownames(summ0)) {
                bhat0[j] <- summ0[ixstr,1]
                pval0[j] <- summ0[ixstr,4]  
                } else {
                bhat0[j] <- pval0[j] <- NA
                }                        
              aic0[j] <- modFit(mod0,selection)
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
              currentAIC[xtest[i]] <- modFit(mod0,selection)
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
    modOK$call <- auxcall
    }
  modOK
  }

# compute lag effects of a covariate (internal use only)
lagEff <- function(model,x=NULL,cumul=F,conf=0.95) {
  #if(class(model)!="lm") stop("Argument 'model' must be an object of class 'lm'",call.=F)
  if(is.null(x)) stop('The name of the covariate is missing',call.=F)
  formstr <- strsplit(strsplit(model$call,"~")[[1]][2],"\\+")[[1]]
  #####
  if(grepl(paste("quec\\(",x,",",sep=""),model$call)) {
    xstr <- formstr[which(grepl(paste("quec\\(",x,",",sep=""),formstr))]
    sx <- as.numeric(strsplit(xstr,",")[[1]][2])
    dx <- as.numeric(strsplit(strsplit(xstr,",")[[1]][3],")")[[1]][1])
    xgrid <- seq(0,dx)
    xtheta <- paste("theta0_quec.",x,sep="")
    imu <- model$coefficients[xtheta]
    icov <- model$vcov[xtheta,xtheta]
    iH <- c()
    for(i in 1:length(xgrid)) {
      if(xgrid[i]<sx-1 | xgrid[i]>dx+1) {
        iH[i] <- 0
        } else {
        iH[i] <- (-4/(sx-dx-2)^2)*(xgrid[i]^2-(sx+dx)*xgrid[i]+(sx-1)*(dx+1))
        }
      }
    iH <- matrix(iH,ncol=1)
    #####
    } else if(grepl(paste("qdec\\(",x,",",sep=""),model$call)) {
    xstr <- formstr[which(grepl(paste("qdec\\(",x,",",sep=""),formstr))]
    sx <- as.numeric(strsplit(xstr,",")[[1]][2])
    dx <- as.numeric(strsplit(strsplit(xstr,",")[[1]][3],")")[[1]][1])
    xgrid <- seq(0,dx)
    xtheta <- paste("theta0_qdec.",x,sep="")
    imu <- model$coefficients[xtheta]
    icov <- model$vcov[xtheta,xtheta]
    iH <- c()
    for(i in 1:length(xgrid)) {
      if(xgrid[i]<sx | xgrid[i]>dx) {
        iH[i] <- 0
        } else {
        iH[i] <- (xgrid[i]^2-2*(dx+1)*xgrid[i]+(dx+1)^2)/(dx-sx-2)^2
        }
      }
    iH <- matrix(iH,ncol=1)
    #} else if(grepl(paste("trap\\(",x,",",sep=""),model$call)) {
    #xstr <- formstr[which(grepl(paste("trap\\(",x,",",sep=""),formstr))]
    #sx <- as.numeric(strsplit(xstr,",")[[1]][2])
    #isxlag <- as.numeric(strsplit(xstr,",")[[1]][3])
    #idxlag <- as.numeric(strsplit(xstr,",")[[1]][4])
    #dx <- as.numeric(strsplit(strsplit(xstr,",")[[1]][5],")")[[1]][1])
    #xgrid <- seq(0,dx)
    #xtheta <- paste("theta0_trap.",x,sep="")
    #imu <- model$coefficients[xtheta]
    #icov <- model$vcov[xtheta,xtheta]
    #iH <- c()
    #for(i in 1:length(xgrid)) {
    #  if(xgrid[i]>sx-1 & xgrid[i]<isxlag) {
    #    iH[i] <- (xgrid[i]-sx+1)/(isxlag-sx+1)
    #    } else if(xgrid[i]>=isxlag & xgrid[i]<idxlag) {
    #    iH[i] <- 1
    #    } else if(xgrid[i]>=idxlag & xgrid[i]<dx+1) {
    #    iH[i] <- (xgrid[i]-dx-1)/(idxlag-dx-1)
    #    } else {
    #    iH[i] <- 0
    #    }
    #  }
    #iH <- matrix(iH,ncol=1)
    } else if(grepl(paste("gamma\\(",x,",",sep=""),model$call)) {
    xstr <- formstr[which(grepl(paste("gamma\\(",x,",",sep=""),formstr))]
    delta <- as.numeric(strsplit(xstr,",")[[1]][2])
    lambda <- as.numeric(strsplit(strsplit(xstr,",")[[1]][3],")")[[1]][1])
    sx <- 0
    mingam <- qgamma(0.01,1/(1-delta),-log(lambda))+sx
    maxgam <- qgamma(0.99,1/(1-delta),-log(lambda))+sx
    xtheta <- paste("theta0_gamma.",x,sep="")
    xgrid <- seq(0,maxgam)
    imu <- model$coefficients[xtheta]
    icov <- model$vcov[xtheta,xtheta]
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
    if(grepl(x,model$call)==F) stop("Unknown covariate ",x,call.=F)
    xgrid <- 0  
    imu <- model$coefficients[x]
    iH <- matrix(1,nrow=1,ncol=1)
    icov <- matrix(diag(rep(1,length(imu))),nrow=length(imu),ncol=length(imu))
    }                                    
  ibhat <- iH%*%imu
  ibse <- sqrt(diag(iH%*%icov%*%t(iH)))
  quan <- -qnorm((1-conf)/2)
  out <- cbind(ibhat-quan*ibse,ibhat,ibhat+quan*ibse)
  if(cumul==F) {
    rownames(out) <- xgrid
    colnames(out) <- c(paste(100*(1-conf)/2,"%",sep=""),"50%",paste(100*(1+conf)/2,"%",sep=""))
    out
    } else {
    outC <- out
    for(j in 1:ncol(out)) {
      for(i in 1:nrow(out)) {
        outC[i,j] <- sum(out[1:i,j])
        }
      }
    rownames(outC) <- xgrid
    colnames(outC) <- c(paste(100*(1-conf)/2,"%",sep=""),"50%",paste(100*(1+conf)/2,"%",sep=""))
    outC
    }
  }

# adf test (internal use only)
adft <- function(x,k=trunc((length(x)-1)^(1/3))) {
  k <- k + 1
  x <- as.vector(x,mode="double")
  y <- diff(x)
  n <- length(y)
  z <- embed(y,k)
  yt <- z[,1]
  xt1 <- x[k:n]
  tt <- k:n
  if (k > 1) {
    yt1 <- z[,2:k]
    res <- lm(yt ~ xt1 + 1 + tt + yt1)
    }
  else res <- lm(yt ~ xt1 + 1 + tt)
  res.sum <- summary(res)
  STAT <- res.sum$coefficients[2,1]/res.sum$coefficients[2,2]
  table <- cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96), c(3.95, 
      3.8, 3.73, 3.69, 3.68, 3.66), c(3.6, 3.5, 3.45, 3.43, 
      3.42, 3.41), c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12), c(1.14, 
      1.19, 1.22, 1.23, 1.24, 1.25), c(0.8, 0.87, 0.9, 0.92, 
      0.93, 0.94), c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66), c(0.15, 
      0.24, 0.28, 0.31, 0.32, 0.33))
  table <- -table
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
  e <- residuals(lm(x ~ t))
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
 
# function to plot a lag shape (internal use only)
makeShape <- function(bmat,tipo,maxlag=NULL,cumul=F,conf=0.95,ylim=NULL,title=NULL) {    
  if(!is.null(maxlag)) {
    if(length(maxlag)!=1 || !is.numeric(maxlag) || maxlag<0 || maxlag!=round(maxlag)) stop("Argument 'maxlag' must be a non-negative integer number",call.=F)
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
    } else {
    if(cumul==F) {
      bmat <- rbind(bmat,rep(0,ncol(bmat)))
      } else {
      bmat <- rbind(bmat,bmat[nrow(bmat),])
      }
    }                          
  xaux <- (1:nrow(bmat))-2
  upLim <- 1.05*max(bmat)
  if(is.null(ylim)) {
    lowLim <- 1.05*min(bmat)
    upLim <- max(abs(c(upLim,lowLim)))
    lowLim <- -max(abs(c(upLim,lowLim)))
    } else {
    if(length(ylim)!=2) stop("Invalid 'ylim'",call.=F)
    lowLim <- ylim[1]
    upLim <- ylim[2]
    }
  plot(0,type="n",xlim=c(xaux[1],rev(xaux)[1]),ylim=c(lowLim,upLim),yaxs="i",xaxs="i",cex.lab=1.2,
    lwd=2,xaxt="n",yaxt="n",xlab="Lag",ylab="Coefficient",main=title,cex.main=1.2)
  if(cumul==T) {
    mtext("cumulative lag shape",cex=0.9)
    } else {
    mtext("instantaneous lag shape",cex=0.9)    
    }
  usesplin <- 0
  auxsplin <- bmat[which(bmat[,2]!=0),2]
  if(length(auxsplin)>=3) usesplin <- 1
  if(usesplin==0) {
    xgrid <- xaux
    ygrid <- bmat
    } else {
    xgrid <- seq(min(xaux[which(bmat[,2]!=0)])-1,max(xaux[which(bmat[,2]!=0)])+1,length=1000)
    ygrid <- matrix(nrow=length(xgrid),ncol=ncol(bmat))
    splind <- max(1,(min(which(bmat[,2]!=0))-1)):min(length(xaux),(max(which(bmat[,2]!=0))+1))
    for(i in 1:ncol(bmat)) {
      auxspl <- smooth.spline(xaux[splind],bmat[splind,i])
      ygrid[,i] <- predict(auxspl,xgrid)$y
      }
    }
  polygon(c(xgrid,rev(xgrid)),c(ygrid[,1],rev(ygrid[,3])),border=NA,col="grey83")
  yaxaux <- seq(lowLim,upLim,length=21)
  ylabaux <- signif(yaxaux,3)
  ylabaux[11] <- 0
  xaxaux <- seq(min(xaux),max(xaux))
  auxby <- max(1,round((max(xaux)-min(xaux)+1)/30))
  xlabaux1 <- xlabaux2 <- seq(min(xaux),max(xaux),by=auxby)
  xlabaux2[c(1,length(xlabaux1))] <- NA
  abline(h=yaxaux,v=seq(min(xaux),max(xaux),by=auxby),col="grey75",lty=2)
  abline(h=0,lty=2,col="grey35")
  lines(ygrid[,2]~xgrid,col="grey35",lty=2)
  auxpoi1 <- which(bmat[,2]!=0)
  auxpoi2 <- max(1,min(auxpoi1)-1):min(length(xaux),max(auxpoi1)+1)
  points(bmat[auxpoi2,2]~xaux[auxpoi2],col="grey35",lty=2,cex=0.5)
  axis(1,at=xlabaux1,labels=xlabaux2,cex.axis=1.1)
  axis(2,at=yaxaux,labels=ylabaux,cex.axis=1.1)       
  if(cumul==F) {
    tcaVal <- signif(apply(bmat,2,sum),5)
    } else {
    tcaVal <- signif(bmat[nrow(bmat),],5)
    }
  confLeg <- paste("   ",conf*100,"% CI: (",tcaVal[1],", ",tcaVal[3],")",sep="")
  if(max(bmat[,2])>0) {
    legpos <- "bottomright"
    } else {
    legpos <- "topright"
    }                                       
  est <- bmat[,2]
  if(cumul==T) {
    newest <- est[1]
    for(i in 2:length(est)) {
      newest[i] <- est[i]-est[i-1]
      } 
    est <- newest
    }
  minlag <- min(as.numeric(rownames(bmat)[which(est!=0)]))
  maxlag <- max(as.numeric(rownames(bmat)[which(est!=0)]))
  legend(legpos,legend=c(paste("Effective lags: ",minlag," to ",maxlag,sep=""),paste("Cumulative coefficient: ",tcaVal[2],sep=""),confLeg),bty="n",cex=1.1)
  box()
  }    

# check if a vector is named (internal use only)
hasname <- function(x) {
  sum(sapply(names(x),function(z){nchar(z)>0}))==length(x)
  }
  
# check control options (internal use only)
controlCheck <- function(parstr,control,parSets,is.sign=F) {
  if(parstr %in% names(control)) {
    if(is.sign==F) {
      if(!is.list(control[[parstr]]) || hasname(control[[parstr]])==F || sum(sapply(control[[parstr]],class)!="numeric")>0 || sum(sapply(control[[parstr]],function(x){x!=round(x)||x<0}))>0 ) {
        stop("Component '",parstr,"' of argument 'control' must be a named list of non-negative integers",call.=F)
        }
      } else {
      if(!is.list(control[[parstr]]) || hasname(control[[parstr]])==F || sum(sapply(control[[parstr]],function(x){sum(x %in% c("+","-"))<length(x)})>0)) {
        stop("Component '",parstr,"' of argument 'control' must be a named list of characters '+' or '-'",call.=F)
        }
      }
    }
  auxch <- setdiff(names(control[[parstr]]),names(parSets))
  if(length(auxch)>0) stop("Unknown variables in component '",parstr,"' of argument 'control': ",paste(auxch,collapse=", "),call.=F)
  auxG <- lapply(control[[parstr]],names)
  gnam <- names(auxG)
  for(i in 1:length(gnam)) {
    auxch <- setdiff(auxG[[gnam[i]]],parSets[[gnam[i]]])
    if(length(auxch)>0) stop("Unknown covariates of variable ",gnam[i]," in component '",parstr,"' of argument 'control': ",paste(auxch,collapse=", "),call.=F)
    }
  
  }

# fit a distributed-lag SEM
dlsem <- function(model.code,group=NULL,exogenous=NULL,data,log=FALSE,control=NULL,imputation=TRUE,uniroot.check=TRUE,test="adf",combine="choi",k=0,lshort=TRUE,maxdiff=5,tol=0.0001,maxit=500,selection="aic",plotDir=NULL) {
  if(sum(sapply(model.code,class)!="formula")>0) stop("Argument 'model code' must be a list of formulas",call.=F)
  if(!is.null(control) && !is.list(control)) stop("Argument 'control' must be a list",call.=F) 
  if(imputation==T) {  
    if(length(maxit)!=1 || maxit<=0 || maxit!=round(maxit)) stop("Argument 'maxit' must be a non-negative integer number",call.=F)
    if(length(tol)!=1 || tol<=0) stop("Argument 'tol' must be a positive real number",call.=F) 
    }
  if(uniroot.check==T) {
    if(length(test)!=1 || (test %in% c("adf","kpss"))==F) stop("Argument 'test' can be either 'adf' or 'kpss'",call.=F)
    if(length(combine)!=1 || (combine %in% c("choi","demetrescu"))==F) stop("Argument 'combine' can be either 'choi' or 'demetrescu'",call.=F)
    if(length(maxdiff)!=1 || maxdiff<0 || maxdiff!=round(maxdiff)) stop("Argument 'maxdiff' must be a non-negative integer number",call.=F)
    }
  if((selection %in% c("aic","bic","mdl"))==F) stop("Argument selection must be one among 'aic', 'bic' and 'mdl'",call.=F)
  res <- list()
  messlen <- 0
  pset <- list()
  for(i in 1:length(model.code)) {
    icheck <- scanForm(model.code[[i]])
    if(length(icheck$lpar)>0) {
      for(j in 1:length(icheck$lpar)) {
        if(icheck$ltype[j]=="quec") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>=icheck$lpar[[j]][2]) {
            stop("Invalid parameter values for lag shape 'quec' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="qdec") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>=icheck$lpar[[j]][2]) {
            stop("Invalid parameter values for lag shape 'qdec' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="gamma") {
          if(length(icheck$lpar[[j]])!=2 || sum(icheck$lpar[[j]]<=0)>0 || sum(icheck$lpar[[j]]>=1)>0) {
            stop("Invalid parameter values for lag shape 'gamma' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          #} else if(icheck$ltype[j]=="trap") {
          #if(length(icheck$lpar[[j]])!=4 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 ||
          #  icheck$lpar[[j]][1]>icheck$lpar[[j]][2] || icheck$lpar[[j]][2]>icheck$lpar[[j]][3] || icheck$lpar[[j]][3]>icheck$lpar[[j]][4] || icheck$lpar[[j]][4]-icheck$lpar[[j]][1]<=0) {
          #  stop("Invalid parameter values for lag shape 'trap' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
          #  }
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
    auxdupl <- duplicated(xnam)
    if(sum(auxdupl)>0) {
      stop("Duplicated covariate for ",ynam,": ",xnam[auxdupl][1],call.=F)
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
  auxvar <- setdiff(c(nodenam,exogenous),colnames(data))
  if(length(auxvar)>0) {
    data[,auxvar] <- as.numeric(rep(NA,nrow(data)))
    cat(paste("Unobserved variables: ",paste(auxvar,collapse=", "),sep=""),"\n")
    imputation <- T
    }
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
  auxch <- setdiff(names(control),c("L","adapt","max.gestation","min.width","max.width","sign"))
  if(length(auxch)>0) stop("Unknown component '",auxch[1],"' in argument 'control'",call.=F)
  if("L" %in% names(control)) {
    if(!is.numeric(control[["L"]]) | hasname(control[["L"]])==F) stop("Component 'L' of argument 'control' must be a named numeric vector",call.=F)   
    auxch <- setdiff(names(control[["L"]]),topG)
    if(length(auxch)>0) stop("Unknown variables in component 'L' of argument 'control': ",paste(auxch,collapse=", "),call.=F)
    }
  if("adapt" %in% names(control)) {
    if(!is.logical(control[["adapt"]]) | hasname(control[["adapt"]])==F) stop("Component 'adapt' of argument 'control' must be a named logical vector",call.=F)
    auxch <- setdiff(names(control[["adapt"]]),topG)
    if(length(auxch)>0) stop("Unknown variables in component 'adapt' of argument 'control': ",paste(auxch,collapse=", "),call.=F)
    }
  if("max.gestation" %in% names(control)) controlCheck("max.gestation",control,pset)
  if("min.width" %in% names(control)) controlCheck("min.width",control,pset)
  if("max.width" %in% names(control)) controlCheck("max.width",control,pset)
  if("sign" %in% names(control)) controlCheck("sign",control,pset,is.sign=T)
  auxch <- intersect(names(control[["min.width"]]),names(control[["max.width"]]))  
  if(length(auxch)>0) {
    for(i in 1:length(auxch)) {
      iminw <- control[["min.width"]][[auxch[i]]]
      imaxw <- control[["max.width"]][[auxch[i]]]
      ich <- intersect(names(iminw),names(imaxw))
      if(length(ich)>0) {
        for(j in 1:length(ich)) {
          if(iminw[ich[j]]>imaxw[ich[j]]) stop("Incoherent components 'min.width' and 'max.width' in argument 'control': ",auxch[i],"~",ich[j],call.=F)
          }
        }
      }
    }
  nodenam <- c(exogenous,topG)
  origdat <- data
  if(log==T) {
    nolog <- c()
    for(i in 1:length(nodenam)) {
      if(sum(data[,nodenam[i]]<=0,na.rm=T)==0) {
        data[,nodenam[i]] <- log(data[,nodenam[i]])
        } else {
        nolog <- c(nolog,nodenam[i])
        }
      }                                                                                             
    if(length(nolog)>0) warning("Logarithmic transformation not applied to variables: ",paste(nolog,collapse=", "),call.=F)
    }
  if(uniroot.check==T) {
    fine <- ndiff <- 0
    cat("Checking stationarity...")
    flush.console()
    while(fine==0) {
      auxp <- c()
      for(i in 1:length(nodenam)) {
        auxdat <- na.omit(data[,nodenam[i]])
        ipvl <- unirootTest(nodenam[i],group=group,data=data,test=test,combine=combine,k=k,lshort=lshort)$p.value
        ifelse(is.null(ipvl), auxp[i] <- 0, auxp[i] <- ipvl)
        }   
      if(sum(auxp>0.05)>0) {
        if(ndiff<maxdiff) {
          data <- applyDiff(nodenam,group=group,data=data,k=1)                                        
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
    }
  if(imputation==T) {                   
    auxna <- apply(data[,nodenam],1,function(x){sum(is.na(x))})
    if(sum(auxna)>0) {
      auxOK <- unname(which(auxna<length(nodenam))) 
      if(uniroot.check==T) {
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
        }
      if(sum(is.na(data[auxOK,]))>0) data[auxOK,] <- EM.imputation(x=nodenam,group=group,data=data[auxOK,],tol=tol,maxit=maxit)
      }
    }
  nomi <- c()
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
    if("L" %in% names(control)) {
      if(nomi[i] %in% names(control[["L"]])) {
        iL <- control[["L"]][nomi[i]]
        } else {
        iL <- 0
        }
      } else {
      iL <- 0
      }
    if("adapt" %in% names(control)) {
      if(nomi[i] %in% names(control[["adapt"]])) {
        iad <- control[["adapt"]][nomi[i]]
        } else {
        iad <- F
        }
      } else {
      iad <- F
      }
    if("max.gestation" %in% names(control)) {
      if(nomi[i] %in% names(control[["max.gestation"]])) {
        iges <- control[["max.gestation"]][[nomi[i]]]
        } else {
        iges <- NULL
        }
      } else {
      iges <- NULL
      }
    if("min.width" %in% names(control)) {
      if(nomi[i] %in% names(control[["min.width"]])) {
        ilen <- control[["min.width"]][[nomi[i]]]
        } else {
        ilen <- NULL
        }
      } else {
      ilen <- NULL
      }
    if("sign" %in% names(control)) {
      if(nomi[i] %in% names(control[["sign"]])) {
        isg <- control[["sign"]][[nomi[i]]]
        } else {
        isg <- NULL
        }
      } else {
      isg <- NULL
      }
    mod0 <- dlaglm(iform,group=group,data=data,log=F,L=iL,adapt=iad,max.gestation=iges,min.width=ilen,sign=isg,selection=selection)
    res[[i]] <- mod0
    messlen <- nchar(auxmess)
    if(!is.null(plotDir)) {
      formstr <- strsplit(strsplit(mod0$call,"~")[[1]][2],"\\+")[[1]]
      wh1 <- which(grepl("quec\\(",formstr))
      wh2 <- which(grepl("qdec\\(",formstr))
      #wh3 <- which(grepl("trap\\(",formstr))
      wh4 <- which(grepl("gamma\\(",formstr))
      ilagged <- c()
      if(length(wh1)>0) {
        for(w in 1:length(wh1)) ilagged <- c(ilagged,strsplit(strsplit(formstr[wh1[w]],"quec\\(")[[1]][2],",")[[1]][1])
        }
      if(length(wh2)>0) {
        for(w in 1:length(wh2)) ilagged <- c(ilagged,strsplit(strsplit(formstr[wh2[w]],"qdec\\(")[[1]][2],",")[[1]][1])
        }
      #if(length(wh3)>0) {
      #  for(w in 1:length(wh3)) ilagged <- c(ilagged,strsplit(strsplit(formstr[wh3[w]],"trap\\(")[[1]][2],",")[[1]][1])
      #  }
      if(length(wh4)>0) {
        for(w in 1:length(wh4)) ilagged <- c(ilagged,strsplit(strsplit(formstr[wh4[w]],"gamma\\(")[[1]][2],",")[[1]][1])
        }
      if(length(ilagged)>0) {
        for(j in 1:length(ilagged)) {
          yaux <- rbind(rep(0,3),lagEff(mod0,ilagged[j],cumul=F,conf=0.95))    
          rownames(yaux)[1] <- "-1"
          if(!is.null(yaux)) {
            pdf(file.path(plotDir,paste(nomi[i],"~",ilagged[j],".pdf",sep="")))
            makeShape(yaux,tipo=0,cumul=F,conf=0.95,title=paste(nomi[i]," ~ ",ilagged[j],sep=""))
            dev.off()
            }
          }
        }
      }
    }
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
    iscan <- scanForm(res[[topG[i]]]$call)
    icov <- setdiff(iscan$X,c(group,exogenous))
    if(length(icov)==0) icov <- "1"
    finalcode[[i]] <- formula(paste(iscan$y," ~ ",paste(icov,collapse="+"),sep=""))
    }
  names(finalcode) <- topG
  out <- list(estimate=res[topG],model.code=finalcode,exogenous=exogenous,group=group,log=log,ndiff=ndiff,data.orig=origdat[,c(group,nodenam)],data.used=data[,c(group,nodenam)])
  class(out) <- "dlsem"
  out
  }

# print method for class dlsem
print.dlsem <- function(x,...) {
  cat("A distributed-lag structural equation model","\n")
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
    cat(" ",n.e,"/",N.e," significant edges at 95% level","\n",sep="")
    } else {
    cat(" No significant edges at 95% level","\n")  
    }
  }

# summary method for class dlsem
summary.dlsem <- function(object,...) {
  lapply(object$estimate,summary)
  }

# fitted method for class dlsem
fitted.dlsem <- function(object,...) {
  pred <- list()
  for(i in 1:length(object$estimate)) {
    pred[[i]] <- object$estimate[[i]]$fitted.values
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
  res
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
  res
  }

# predict method for class dlsem
predict.dlsem <- function(object,...) {
  pred <- lapply(object$estimate,predict,...)
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
  ##### check argument 'newdata'
  res
  }

# plot the lag shape associated to an overall causal effect or a path
lagPlot <- function(model,from=NULL,to=NULL,path=NULL,maxlag=NULL,cumul=FALSE,conf=0.95,ylim=NULL,title=NULL) {
if(class(model)!="dlsem") stop("Argument 'model' must be an object of class 'dlsem'",call.=F)     
  if(is.null(path)) {                                         
    if(is.null(from)) stop("Argument 'from' is missing",call.=F)
    if(is.null(to)) stop("Argument 'to' is missing",call.=F)
    auxpa <- causalEff(model,from=from,to=to,lag=NULL,cumul=cumul,conf=conf)$overall  
    if(!is.null(auxpa)) {
      auxind <- max(which(auxpa[,2]!=0))
      auxpa <- auxpa[1:auxind,]
      }
    } else {
    auxstr <- strsplit(path,"\\*")[[1]]
    if(!is.null(from) && !identical(from,auxstr[1])) stop("Arguments 'from' and 'path' are in conflict",call.=F)
    if(!is.null(to) && !identical(to,rev(auxstr)[1])) stop("Arguments 'to' and 'path' are in conflict",call.=F)
    from <- auxstr[1]
    to <- rev(auxstr)[1]                                                                                              
    auxpa <- causalEff(model,from=from,to=to,cumul=cumul,conf=conf)       
    if((path %in% names(auxpa))==F) stop("Non significant path. Try to reduce the confidence level",call.=F)
    auxpa <- auxpa[[path]]  
    }                                                
  if(!is.null(auxpa)) {
    yaux <- rbind(rep(0,3),auxpa)
    rownames(yaux)[1] <- "-1"
    } else {
    yaux <- NULL
    }                
  if(is.null(title)) {
    if(is.null(path)) {
      title <- paste(to," ~ ",from,sep="")
      } else {
      title <- paste(auxstr,collapse=" * ")
      }
    }
  if(is.null(yaux)) {
    stop("No significant paths found between the selected variables. Try to reduce the confidence level",call.=F)
    } else {
    makeShape(bmat=yaux,tipo=1,maxlag=maxlag,cumul=cumul,conf=conf,ylim=ylim,title=title)
    }
  }

# differentiation
applyDiff <- function(x=NULL,group=NULL,time=NULL,data,k=0) {
  if(class(data)!="data.frame") stop("Argument 'data' must be a data.frame",call.=F)
  auxvar <- setdiff(c(x,group,time),colnames(data))
  if(length(auxvar)>0) stop(paste("Unknown variable: ",paste(auxvar,collapse=", "),sep=""),call.=F)
  if(is.null(x)) x <- setdiff(colnames(data),c(group,time))
  if(!is.null(time) && !is.numeric(data[,time])) stop("'time' must be a numerical variable",call.=F)  
  if(length(k)==1) {
    if(k<0 | k!=round(k)) stop("Invalid difference order",call.=F)
    k <- rep(k,length(x))
    } else {
    if(length(k)!=length(x) || sum(k<0)>0 | sum(k!=round(k))>0) stop("Invalid difference order",call.=F) 
    }                            
  if(!is.null(time)) data <- data[order(data[,time]),]
  #####
  deltaFun <- function(z,k) {
    if(k>0) {
      auxz <- z
      z[1:k] <- NA
      for(i in (k+1):length(z)) {
        z[i] <- auxz[i]-auxz[i-k]
        }
      }
    z
    }
  #####
  diffdat <- data
  if(is.null(group)) {
    for(w in 1:length(x)) {
      if(is.numeric(data[,x[w]])) {
        if(k[w]>=nrow(data)) stop("Too large difference order for variable ",x[w],call.=F)
        diffdat[,x[w]] <- deltaFun(data[,x[w]],k[w])
        }
      }
    } else {
    data[,group] <- factor(data[,group])
    gruppi <- unique(data[,group])
    for(i in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[i])
      for(w in 1:length(x)) {
        if(is.numeric(data[,x[w]])) {
          if(k[w]>=length(auxind)) stop("Too large difference order for variable ",x[w],call.=F)
          diffdat[auxind,x[w]] <- deltaFun(data[auxind,x[w]],k[w])
          }
        }
      } 
    }
  if(!is.null(time) & !is.null(group)) diffdat <- diffdat[order(data[,group]),]
  diffdat
  }

# ADF panel test
unirootTest <- function(x,group=NULL,time=NULL,data,test="adf",combine="choi",k=0,lshort=TRUE) {
  if(length(test)!=1 || (test %in% c("adf","kpss"))==F) stop("Argument 'test' can be either 'adf' or 'kpss'",call.=F)
  if(length(combine)!=1 || (combine %in% c("choi","demetrescu"))==F) stop("Argument 'combine' can be either 'choi' or 'demetrescu'",call.=F)
  if(class(data)!="data.frame") stop("Argument 'data' must be a data.frame",call.=F)
  auxvar <- setdiff(c(x,group,time),colnames(data))
  if(length(auxvar)>0) stop(paste("Unknown variable: ",paste(auxvar,collapse=", "),sep=""),call.=F)
  if(length(x)!=1) stop("Argument 'x' must contain a signle variable name",call.=F)
  if(!is.null(time) && !is.numeric(data[,time])) stop("'time' must be a numerical variable",call.=F)
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

# missing imputation using EM
EM.imputation <- function(x=NULL,group=NULL,data,plotDir=NULL,tol=0.0001,maxit=500,quiet=FALSE) {
  if(class(data)!="data.frame") stop("Argument 'data' must be a data.frame",call.=F)
  if(is.null(colnames(data))) stop("Data must have column names",call.=F)
  if(length(tol)!=1 || tol<=0) stop("Argument 'tol' must be a positive real number",call.=F) 
  if(length(maxit)!=1 || maxit<=0 || maxit!=round(maxit)) stop("Argument 'maxit' must be a non-negative integer number",call.=F)
  if(is.null(x)) x <- setdiff(colnames(data),group)
  auxcheck <- setdiff(c(x,group),colnames(data))
  if(length(auxcheck)>0) stop("Unknown variable: ",paste(auxcheck,collapse=", "),call.=F)           
  for(i in 1:length(x)) {
    if(!is.numeric(data[,x[i]])) stop("Variable ",x[i]," is not numeric",call.=F)
    }
  if(length(x)>1) {                     
    startMu <- apply(data[,x],2,mean,na.rm=T)
    isna <- which(is.nan(startMu))
    if(length(isna)>0) startMu[isna] <- 0
    startS <- matrix(nrow=length(x),ncol=length(x))
    for(i in 1:length(x)) {
      for(j in i:length(x)) {
        auxcov <- na.omit(data[,c(x[i],x[j])]) 
        if(nrow(auxcov)>0) {
          startS[i,j] <- startS[j,i] <- cov(auxcov)[1,2]
          } else {
          startS[i,j] <- startS[j,i] <- 1  #####         
          }
        }    
      }
    startS <- as.matrix(nearPD(startS)$mat)
    if(is.null(group)) {
      old.mu <- startMu
      names(old.mu) <- x
      } else {
      data[,group] <- factor(data[,group])
      gruppi <- unique(data[,group])
      old.mu <- matrix(NA,nrow=length(gruppi),ncol=length(x))  
      for(i in 1:nrow(old.mu)) {
        old.mu[i,] <- startMu
        }
      rownames(old.mu) <- gruppi
      colnames(old.mu) <- x 
      }
    old.S <- startS              
    colnames(old.S) <- rownames(old.S) <- x             
    old.lik <- -Inf
    fine <- 0
    count <- 1                            
    if(quiet==F) cat("Starting EM...")
    flush.console()
    while(fine==0) {
      auxdat <- data
      for(i in 1:nrow(data)) {
        for(j in 1:length(x)) {
          if(is.na(data[i,x[j]])) {                 
            xaux <- x[which(!is.na(data[i,x]))]                   
            if(is.null(group)) {
              jmu <- old.mu
              } else {
              jmu <- old.mu[as.character(data[i,group]),]
              }
            if(length(xaux)>0) {                                      
              xval <- unlist(data[i,xaux]) 
              auxdat[i,x[j]] <- c(jmu[x[j]]+old.S[x[j],xaux]%*%solve(old.S[xaux,xaux])%*%(xval-jmu[xaux]))   
              } else {                    
              auxdat[i,x[j]] <- jmu[x[j]]
              }                                       
            }                 
          }                        
        }
      if(is.null(group)) {
        new.mu <- colMeans(auxdat[,x])
        } else {                        
        new.mu <- matrix(nrow=length(gruppi),ncol=length(x))  
        rownames(new.mu) <- gruppi
        colnames(new.mu) <- x 
        for(i in 1:length(gruppi)) {                          
          new.mu[i,] <- colMeans(auxdat[which(data[,group]==gruppi[i]),x])
          }
        } 
      new.S <- as.matrix(nearPD(cov(auxdat[,x]))$mat)
      new.lik <- 0                           
      iS <- solve(new.S)
      for(i in 1:nrow(data)) {              
        if(is.null(group)) {
          idel <- unlist(auxdat[i,x])-new.mu
          } else {
          idel <- unlist(auxdat[i,x])-new.mu[as.character(data[i,group]),]           
          }                                                                                  
        new.lik <- new.lik-0.5*length(x)*log(2*pi)-0.5*log(det(new.S))-0.5*t(idel)%*%iS%*%idel
        }
      auxmess <- paste("EM iteration ",count,". Log-likelihood: ",round(new.lik,4),sep="")
      if(quiet==F) cat('\r',auxmess,sep="")       
      flush.console()                                   
      if(new.lik-old.lik<tol | count>maxit) {
        fine <- 1
        if(count<maxit) {
          auxfinal <- paste("EM converged after ",count-1," iterations. Log-likelihood: ",round(new.lik,4),sep="")
          } else {
          auxfinal <- paste("Maximum number of iterations exceeded. Log-likelihood: ",round(new.lik,4),sep="")      
          }
        nchardiff <- nchar(auxmess)-nchar(auxfinal)
        if(nchardiff>0) {
          delta <- rep(" ",nchardiff)
          } else {
          delta <- ""
          }
        if(quiet==F) cat('\r',auxfinal,delta,sep="","\n")
        } else {
        old.lik <- new.lik
        old.mu <- new.mu
        old.S <- new.S
        count <- count+1
        }
      }                                         
    } else {
    auxdat <- data
    if(is.null(group)) {
      if(sum(!is.na(data[,x]))>=2) auxdat[which(is.na(data[,x])),x] <- mean(data[,x],na.rm=T)
      } else {
      data[,group] <- factor(data[,group])
      gruppi <- unique(data[,group])
      for(i in 1:length(gruppi)) {
        coudat <- data[which(data[,group]==gruppi[i]),x]
        if(sum(!is.na(coudat))>=2) {
          coudat[which(is.na(coudat))] <- mean(coudat,na.rm=T)
          auxdat[which(data[,group]==gruppi[i]),x] <- coudat
          }
        }
      }
    }
  if(!is.null(plotDir)) {
    if(!is.null(group)) {
      gruppi <- unique(data[,group])
      for(j in 1:length(x)) {
        for(w in 1:length(gruppi)) {
          auxind <- which(auxdat[,group]==gruppi[w])
          xvaljw <- cbind(auxdat[auxind,x[j]],data[auxind,x[j]],1:length(auxind))                                 
          timelab <- "ID"
          if(sum(!is.na(xvaljw[,1])>0)) {
            pdf(file.path(plotDir,paste(x[j],"_(",gruppi[w],").pdf",sep="")))
            plot(xvaljw[,1]~xvaljw[,3],xlab=timelab,ylab=x[j],type="n",main=paste(x[j]," (",gruppi[w],")",sep=""),cex.main=1.4) 
            grid(col="grey85",lty=2)
            points(xvaljw[,1]~xvaljw[,3],cex=0.8)
            lines(xvaljw[,1]~xvaljw[,3],cex=0.8,lty=2)
            auxna <- which(is.na(xvaljw[,2]))
            if(length(auxna)>0) points(xvaljw[auxna,1]~xvaljw[auxna,3],cex=0.8,col="red")
            box()
            dev.off()
            }
          }
        }
      } else {
      for(j in 1:length(x)) {
        xvalj <- cbind(auxdat[,x[j]],data[,x[j]],1:nrow(data))                                 
        timelab <- "ID"
        if(sum(!is.na(xvalj[,1])>0)) {
          pdf(file.path(plotDir,paste(x[j],".pdf",sep="")))
          plot(xvalj[,1]~xvalj[,3],xlab=timelab,ylab=x[j],type="n",main=paste(x[j],sep=""),cex.main=1.4) 
          grid(col="grey80",lty=2)
          points(xvalj[,1]~xvalj[,3],cex=0.8)
          lines(xvalj[,1]~xvalj[,3],cex=0.8,lty=2)
          auxna <- which(is.na(xvalj[,2]))
          if(length(auxna)>0) points(xvalj[auxna,1]~xvalj[auxna,3],cex=0.8,col="red")
          box()
          dev.off()
          }
        }
      }
    }
  auxdat  
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
      ijnam <- auxnam[j]
      ijnam <- gsub("theta0_quec.","",ijnam)
      ijnam <- gsub("theta0_qdec.","",ijnam)
      #ijnam <- gsub("theta0_trap.","",ijnam)
      ijnam <- gsub("theta0_gamma.","",ijnam)
      if(ijnam %in% nomi) {        
        G0 <- addEdge(ijnam,nomi[i],G0,1)
        if(isumm[auxnam[j],4]<1-conf) {
          G <- addEdge(ijnam,nomi[i],G,1)
          if(auxsg>0) {
            eSign <- c(eSign,"+")
            } else {
            eSign <- c(eSign,"-")
            }
          } else {
          eSign <- c(eSign,"")
          }
        names(eSign)[length(eSign)] <- paste(ijnam,"~",nomi[i],sep="")              
        }
      }
    }
  list(graph=G,full.graph=G0,sign=eSign)
  }

# compute the coefficient associated to each edge at different time lags
edgeCoeff <- function(x,lag=NULL,conf=0.95) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  nomi <- names(x$estimate)
  if(is.null(lag)) {
    laglen <- c()
    for(i in 1:length(nomi)) {
      isumm <- summary(x$estimate[[nomi[i]]])$coefficients
      auxnam <- rownames(isumm)
      for(j in 1:length(auxnam)) {
        ijnam <- auxnam[j]
        ijnam <- gsub("theta0_quec.","",ijnam)
        ijnam <- gsub("theta0_qdec.","",ijnam)
        #ijnam <- gsub("theta0_trap.","",ijnam)
        ijnam <- gsub("theta0_gamma.","",ijnam)
        if(ijnam %in% nomi) {
          cumb <- lagEff(x$estimate[[nomi[i]]],x=ijnam,cumul=F,conf=conf)
          laglen <- c(laglen,rownames(cumb)[nrow(cumb)])
          }
        }
      }                                         
    lag <- 0:max(as.numeric(laglen))
    }
  bList <- vector("list",length=length(lag))
  for(i in 1:length(nomi)) {
    isumm <- summary(x$estimate[[nomi[i]]])$coefficients
    auxnam <- rownames(isumm)
    for(j in 1:length(auxnam)) {
      ijnam <- auxnam[j]
      ijnam <- gsub("theta0_quec.","",ijnam)
      ijnam <- gsub("theta0_qdec.","",ijnam)
      #ijnam <- gsub("theta0_trap.","",ijnam)
      ijnam <- gsub("theta0_gamma.","",ijnam)
      if(ijnam %in% nomi) {
        cumb <- lagEff(x$estimate[[nomi[i]]],x=ijnam,cumul=F,conf=conf)
        for(w in 1:length(lag)) {
          if(as.character(lag[w]) %in% rownames(cumb)) {
            wbhat <- cumb[as.character(lag[w]),]
            } else {
            wbhat <- rep(0,3)
            }
          if(isumm[auxnam[j],4]<1-conf) {
            bList[[w]] <- rbind(bList[[w]],wbhat)                       
            } else {
            bList[[w]] <- rbind(bList[[w]],rep(0,3))               
            }
          rownames(bList[[w]])[nrow(bList[[w]])] <- paste(ijnam,"~",nomi[i],sep="")
          }
        }
      }
    if(!is.null(bList)) names(bList) <- lag
    }
  for(i in 1:length(bList)) {
    auxnam <- rownames(bList[[i]])
    newnam <- c()
    for(j in 1:length(auxnam)) {
      newnam[j] <- paste(rev(strsplit(auxnam[j],"~")[[1]]),collapse="~")
      }
    rownames(bList[[i]]) <- newnam
    }
  for(i in 1:length(bList)) colnames(bList[[i]]) <- c(paste(100*(1-conf)/2,"%",sep=""),"50%",paste(100*(1+conf)/2,"%",sep=""))
  bList
  }

# plot method for class dlsem
plot.dlsem <- function(x,conf=0.95,node.col=NULL,font.col=NULL,border.col=NULL,edge.col=NULL,edge.lab=NULL,...) {
  defpar <- par()
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
  nAttr$label <- sapply(nomi,cutString,l=12)  ##### max number of characters: 12
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
  #if(sign.col==T) {
  eCol <- G$sign
  eCol[which(G$sign=="+")] <- "green4"
  eCol[which(G$sign=="-")] <- "tomato3"
  eCol[which(G$sign=="")] <- "grey55" #####
  eAttr[[length(eAttr)+1]] <- eCol
  names(eAttr)[length(eAttr)] <- "color"
    #} else {
    #if(!is.null(edge.col)) {
    #  eAttr$color <- edge.col     
    #  eAttr$color[which(G$sign=="")] <- NA
    #  } else {
    #  eCol <- rep("grey25",length(G$sign))
    #  names(eCol) <- names(G$sign)
    #  eCol[which(G$sign=="")] <- NA
    #  eAttr[[length(eAttr)+1]] <- eCol
    #  names(eAttr)[length(eAttr)] <- "color"
    #  }
    #}                              
  plot(G$full.graph,"dot",nodeAttrs=nAttr,edgeAttrs=eAttr,attrs=list(edge=list(color="grey25",arrowsize=0.4)),...)
  options(warn=-1)
  par(defpar)
  options(warn=0)
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
isIndep <- function(x,var1,var2,given=NULL,conf=0.95) {
  if(class(x)!="dlsem") stop("The first argument must be an object of class 'dlsem'",call.=F)
  if(length(var1)!=1) stop("Argument 'var1' must have length 1",call.=F)
  if(length(var2)!=1) stop("Argument 'var2' must have length 1",call.=F)
  G <- makeGraph(x,conf=conf)$graph
  nomi <- nodes(G)
  auxcheck <- setdiff(c(var1,var2,given),nomi)
  if(length(auxcheck)>0) stop("Unknown variable ",paste(auxcheck,collapse=", "),call.=F)
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
  auxlag <- 0:lag
  out <- list()
  for(w in 1:length(auxlag)) {
    mycode <- "res <- c(); "
    for(i in 1:g) {
      mycode <- paste(mycode,"for(k",i," in 1:length(x[[",i,"]])) { ; ",sep="")
      }
    mycode <- paste(mycode,paste("xaux <- c(",paste("x[[",1:g,"]][k",1:g,"]",collapse=",",sep=""),"); ",sep=""),sep="")
    mycode <- paste(mycode,"if(sum(xaux)==auxlag[w]) { res <- rbind(res,xaux) }; ",sep="")
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
  names(out) <- auxlag
  out
  }

# disentanglement of causal effects
causalEff <- function(x,from=NULL,to=NULL,lag=NULL,cumul=FALSE,conf=0.95) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  auxcheck <- setdiff(c(from,to),names(x$estimate))
  if(length(auxcheck)>0) {
    auxcntx <- intersect(auxcheck,x$exogenous)
    if(length(auxcntx)>0) {
      stop("Exogenous variables cannot be involved in path analysis",call.=F)
      } else {
      stop("Unknown variable: ",paste(auxcheck,collapse=", "),call.=F)
      }
    }
  if(is.null(from)) stop("Argument 'from' is missing",call.=F)
  if(is.null(to)) stop("Argument 'to' is missing",call.=F)
  if(!is.character(from)) stop("Invalid argument 'from'",call.=F)
  if(!is.character(to) || length(to)!=1) stop("Invalid argument 'to'",call.=F)
  Gobj <- makeGraph(x,conf=conf)
  G <- Gobj$graph
  nomi <- nodes(G)
  if(length(setdiff(from,nomi))>0) stop("Unknown variable ",paste(setdiff(from,nomi),collapse=", "),call.=F)
  if((to %in% nomi)==F) stop("Unknown variable ",to,call.=F)
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
    if(is.null(lag)) {
      laglen <- list()
      for(i in 1:length(newPathList)) {
        jlaglen <- list()
        for(j in 2:length(newPathList[[i]])) {
          auxnam <- paste(newPathList[[i]][j],"~",newPathList[[i]][j-1],sep="")
          auxeff <- lagEff(x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf)
          auxpos <- which(auxeff[,2]!=0)
          if(length(auxpos)>0) {
            jlaglen[[j-1]] <- as.numeric(rownames(auxeff)[auxpos])
            } else {
            jlaglen[[j-1]] <- NA
            }
          }
        laglen[[i]] <- c(min(jlaglen[[1]],na.rm=T),sum(sapply(jlaglen,max,na.rm=T)))
        }          
      lag <- 0:(max(unlist(laglen))+1)
      cutab <- 0
      } else {
      for(i in 1:length(lag)) {
        if(!is.numeric(lag[i]) || lag[i]<0 || lag[i]!=round(lag[i])) stop("Argument 'lag' must contain non-negative integer numbers only",call.=F)
        }                                       
      lagOK <- lag
      lag <- 0:max(lag)
      cutab <- 1
      }
    mycoeff <- edgeCoeff(x,lag=lag,conf=conf)
    #
    sd_calc <- function(muvet,sdvet) { sqrt(prod(muvet^2+sdvet^2)-prod(muvet^2)) }
    quan <- -qnorm((1-conf)/2)
    bhat <- list()
    for(i in 1:length(mycoeff)) {
      bhat[[i]] <- matrix(nrow=nrow(mycoeff[[i]]),ncol=2)
      for(j in 1:nrow(mycoeff[[i]])) {
        auxeval <- mycoeff[[i]][j,2]
        auxsd <- (mycoeff[[i]][j,3]-mycoeff[[i]][j,2])/quan
        bhat[[i]][j,] <- c(auxeval,auxsd)
        }
      rownames(bhat[[i]]) <- rownames(mycoeff[[i]])
      }
    names(bhat) <- names(mycoeff)
    outList <- list()                       
    for(i in 1:length(newPathList)) {           
      outList[[i]] <- matrix(nrow=length(lag),ncol=3)
      rownames(outList[[i]]) <- lag
      colnames(outList[[i]]) <- c(paste(100*(1-conf)/2,"%",sep=""),"50%",paste(100*(1+conf)/2,"%",sep=""))
      auxbetalag <- list()
      for(j in 2:length(newPathList[[i]])) {
        auxnam <- paste(newPathList[[i]][j],"~",newPathList[[i]][j-1],sep="")
        auxeff <- lagEff(x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf)
        auxpos <- which(auxeff[,2]!=0)
        if(length(auxpos)>0) {                                 
          auxbetalag[[j-1]] <- as.numeric(rownames(auxeff)[auxpos])
          } else {
          auxbetalag[[j-1]] <- 0
          }
        names(auxbetalag)[j-1] <- auxnam
        }                              
      lagsumMat <- findlag2sum(auxbetalag,lag=max(lag))
      for(j in 1:length(lagsumMat)) {
        auxlag <- as.character(lag[j])                           
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
          outList[[i]][j,] <- c(mupath-quan*sdpath,mupath,mupath+quan*sdpath)
          } else {
          outList[[i]][j,] <- rep(0,3)        
          }
        }
      }
    names(outList) <- sapply(newPathList,function(x){paste(x,collapse="*")})
    out <- matrix(nrow=length(lag),ncol=3)
    rownames(out) <- lag
    colnames(out) <- c(paste(100*(1-conf)/2,"%",sep=""),"50%",paste(100*(1+conf)/2,"%",sep=""))
    for(j in 1:length(lag)) {
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
    if(cutab==1) {
      for(i in 1:length(outList)) {
        outList[[i]] <- outList[[i]][as.character(lagOK),]    
        }
      }
    outList
    } else {
    stop("No significant paths found between the selected variables. Try to reduce the confidence level",call.=F)
    }
  }

# fit indices
fitIndices <- function(x) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  xfit <- x$estimate     
  deg <- n <- expdev <- resdev <- Rsq <- aic <- bic <- mdl <- c()
  for(i in 1:length(xfit)) {
    resdev[i] <- deviance(xfit[[i]])
    auxss <- anova(xfit[[i]])$'Sum Sq'
    expdev[i] <- sum(auxss[1:(length(auxss)-1)])
    n[i] <- length(residuals(xfit[[i]]))
    deg[i] <- xfit[[i]]$df.residual
    Rsq[i] <- summary(xfit[[i]])$'r.squared'
    aic[i] <- modFit(xfit[[i]],"aic")
    bic[i] <- modFit(xfit[[i]],"bic")
    mdl[i] <- modFit(xfit[[i]],"mdl")
    }  
  out1 <- c(Rsq,sum(expdev)/(sum(expdev)+sum(resdev)))
  names(out1) <- c(names(xfit),"(overall)")
  out2 <- c(aic,sum(aic))
  names(out2) <- c(names(xfit),"(overall)")
  out3 <- c(bic,sum(bic))
  names(out3) <- c(names(xfit),"(overall)")
  out4 <- c(mdl,sum(mdl))
  names(out4) <- c(names(xfit),"(overall)")
  list('R squared'=out1,'Akaike Information Criterion'=out2,'Bayesian Information Criterion'=out3,'Minimum Description Length'=out4)
  }