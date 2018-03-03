defpar <- par()[setdiff(names(par()),c("cin","cra","csi","cxy","din","page"))]

# lag weights (internal use only)
lagwei <- function(theta,lag,type) {
  res <- c()
  xa <- 0  #####  
  for(i in 1:length(lag)) {
    if(type=="quec.lag") {
      if(lag[i]<theta[1]-1 | lag[i]>theta[2]+1) {
        res[i] <- 0
        } else {
        res[i] <- -4/(theta[2]-theta[1]+2)^2*(lag[i]^2-(theta[1]+theta[2])*lag[i]+(theta[1]-1)*(theta[2]+1))
        }
      } else if(type=="qdec.lag") {
      if(lag[i]<theta[1] | lag[i]>theta[2]+1) {
        res[i] <- 0
        } else {
        res[i] <- (lag[i]^2-2*(theta[2]+1)*lag[i]+(theta[2]+1)^2)/(theta[2]-theta[1]+1)^2
        }
      } else if(type=="gamma.lag") {
      if(lag[i]>=xa) {
        bnum <- (lag[i]-xa+1)^(theta[1]/(1-theta[1]))*theta[2]^(lag[i]-xa)
        xM <- (theta[1]/(theta[1]-1))/log(theta[2])+xa-1
        bden <- (xM-xa+1)^(theta[1]/(1-theta[1]))*theta[2]^(xM-xa)
        res[i] <- bnum/bden
        } else {
        res[i] <- 0     
        }
      }
    }
  names(res) <- lag
  res
  }

# generate lagged instances of a variable (internal use only)
genLag <- function(x,maxlag) {
  if(maxlag>0) {
    out <- x
    for(w in 1:maxlag) {
      if(w<length(x)) {
        wx <- c(rep(NA,w),x[1:(length(x)-w)])
        } else {
        wx <- rep(NA,length(x))
        }
      out <- cbind(out,wx)        
      }
    colnames(out) <- NULL
    out
    } else {
    x
    }
  }

# distributed-lag transformation (internal use only)
Zmat <- function(x,type,theta) {
  if(type=="none") {
    matrix(x,ncol=1)
    } else {
    if(type=="gamma.lag") {
      xa <- 0  #####
      xb <- qgamma(0.99,1/(1-theta[1]),-log(theta[2]))+xa
      laglim <- xb
      } else {
      laglim <- theta[2]
      }
    H <- lagwei(theta,0:laglim,type)
    as.numeric(genLag(x,laglim)%*%matrix(H,ncol=1))
    }
  }

# quec lag transformation (internal use only)
quec.lag <- function(x,a,b,x.group) {
  if(is.null(x.group)) {
    Zmat(x,"quec.lag",c(a,b))
    } else {
    res <- c()  
    gruppi <- levels(factor(x.group))
    for(i in 1:length(gruppi)) {
      auxind <- which(x.group==gruppi[i])
      res[auxind] <- Zmat(x[auxind],"quec.lag",c(a,b))
      }
    res
    }
  }

# qdec lag transformation (internal use only)
qdec.lag <- function(x,a,b,x.group) {
  if(is.null(x.group)) {
    Zmat(x,"qdec.lag",c(a,b))
    } else {
    res <- c()  
    gruppi <- levels(factor(x.group))
    for(i in 1:length(gruppi)) {
      auxind <- which(x.group==gruppi[i])
      res[auxind] <- Zmat(x[auxind],"qdec.lag",c(a,b))
      }
    res
    }
  }

# gamma lag transformation (internal use only)
gamma.lag <- function(x,a,b,x.group) {
  if(is.null(x.group)) {
    Zmat(x,"gamma.lag",c(a,b))
    } else {
    res <- c()  
    gruppi <- levels(factor(x.group))
    for(i in 1:length(gruppi)) {
      auxind <- which(x.group==gruppi[i])
      res[auxind] <- Zmat(x[auxind],"gamma.lag",c(a,b))
      }
    res
    }
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
    gruppi <- levels(factor(data[,group]))
    for(i in 1:length(gruppi)) {
      auxlaglim[i] <- nrow(na.omit(data[which(data[,group]==gruppi[i]),]))
      }
    min(auxlaglim,na.rm=T)
    }
  }

# scan formula (internal use only)
scanForm <- function(x) {
  auxform <- gsub(" ","",formula(x))[-1]
  ynam <- gsub(" ","",auxform[1])
  auX <- gsub(" ","",strsplit(auxform[2],"\\+")[[1]])
  if(sum(!is.na(auX)==0)) auX <- "1"
  auX <- setdiff(auX,c("1","-1"))
  if(length(auX)>0) {
    lnames <- ltype <- rep(NA,length(auX))
    lpar <- list()
    for(i in 1:length(auX)) {        
      if(nchar(auX[i])>0) {
        if(identical("quec.lag(",substr(auX[i],1,9))) {                           
          istr <- gsub("\\)","",gsub("quec\\.lag\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "quec.lag"
          lpar[[i]] <- as.numeric(istr[2:3])
          } else if(identical("qdec.lag(",substr(auX[i],1,9))) {                            
          istr <- gsub("\\)","",gsub("qdec\\.lag\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "qdec.lag"
          lpar[[i]] <- as.numeric(istr[2:3])             
          } else if(identical("gamma.lag(",substr(auX[i],1,10))) {
          istr <- gsub("\\)","",gsub("gamma\\.lag\\(","",strsplit(auX[i],",")[[1]]))
          lnames[i] <- istr[1]
          ltype[i] <- "gamma.lag"
          lpar[[i]] <- as.numeric(istr[2:3])             
          } else {
          lnames[i] <- auX[i]
          ltype[i] <- "none"
          lpar[[i]] <- NA
          }                        
        if(grepl("\\(",lnames[i])) {
          stop("Unknown lag shape: ",strsplit(lnames[i],"\\(")[[1]][1],". Available lag shapes are 'quec.lag', 'qdec.lag' and 'gamma.lag'.",call.=F)
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
  if(lag.type=="gamma.lag") {
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

# create lm formula (internal use only)
creatForm <- function(y,X,group,type,theta) {
  xnam <- c()
  if(length(X)>0) {
    for(i in 1:length(X)) {
      if(type[i]=="none") {
        xnam[i] <- X[i]
        } else {
        if(is.null(group)) {
          xnam[i] <- paste(type[i],"(",X[i],",",theta[[i]][1],",",theta[[i]][2],")",sep="")
          } else {
          xnam[i] <- paste(type[i],"(",X[i],",",theta[[i]][1],",",theta[[i]][2],",",group,")",sep="")
          }
        }
      }
    }
  if(is.null(group)) {
    if(length(X)>0) {
      res <- paste(y,"~",paste(xnam,collapse="+"),sep="")    
      } else {
      res <- paste(y,"~1",sep="")
      }
    } else {
    if(length(X)>0) {
      res <- paste(y,"~-1+",group,"+",paste(xnam,collapse="+"),sep="")    
      } else {
      res <- paste(y,"~-1+",group,sep="")
      }
    }
  formula(res)
  }

# extract the name from a lag shape (inrternal use only)
extrName <- function(x) {
  if(identical("quec.lag(",substr(x,1,9))) {                           
    gsub("\\)","",gsub("quec\\.lag\\(","",strsplit(x,",")[[1]]))[1]
    } else if(identical("qdec.lag(",substr(x,1,9))) {                            
    gsub("\\)","",gsub("qdec\\.lag\\(","",strsplit(x,",")[[1]]))[1]
    } else if(identical("gamma.lag(",substr(x,1,10))) {
    gsub("\\)","",gsub("gamma\\.lag\\(","",strsplit(x,",")[[1]]))[1]
    } else {
    x
    }
  }

# extract lag parameters from a formula (internal use only)
extrLPar <- function(x) {
  auxform <- scanForm(x)
  nomi <- auxform$X
  if(length(nomi)>0) {
    res <- matrix(nrow=length(nomi),ncol=2)
    for(i in 1:length(nomi)) {
      inam <- nomi[i]
      ilab <- extrName(inam)
      if(ilab %in% names(auxform$ltype)) {
        ityp <- auxform$ltype[ilab]
        ipar <- auxform$lpar[[ilab]]
        } else {
        ityp <- "none"
        ipar <- NA
        }
      if(ityp=="none") {
        res[i,] <- c(0,0)
        } else if(ityp=="gamma.lag") {
        idel <- as.character(ipar[1])
        ilam <- as.character(ipar[2])
        igampar <- c(gammaCons$sx[idel,ilam],gammaCons$dx[idel,ilam])
        res[i,] <- c(igampar[1],igampar[2])
        } else {
        res[i,] <- c(ipar[1],ipar[2]) 
        }
      }
    rownames(res) <- sapply(nomi,extrName)
    colnames(res) <- c("a","b") 
    res
    }
  }

# check if qualitative variables have at least 2 levels (internal use only)
checkQV <- function(x,z,group,k,data) {
  dat0 <- data[,c(group,x,z),drop=F]
  x2del <- c()
  if(is.null(group)) {
    auxcheck <- dat0
    auxcheck[1:k,x] <- NA
    auxcheck <- na.omit(auxcheck)
    for(j in 1:length(x)) {
      if(length(unique(auxcheck[,x[j]]))<2) x2del <- c(x2del,x[j])
      }
    } else {
    gruppi <- levels(factor(data[,group]))
    #
    auxind <- which(dat0[,group]==gruppi[1])
    auxcheck <- dat0[auxind,]
    auxcheck[1:k,x] <- NA
    auxcheck <- na.omit(auxcheck)
    for(j in 1:length(x)) {
      if(length(unique(auxcheck[,x[j]]))<2) x2del <- c(x2del,x[j])
      }
    #
    #end <- 0
    #for(i in 1:length(gruppi)) {
    #  auxind <- which(dat0[,group]==gruppi[i])
    #  auxcheck <- dat0[auxind,]
    #  auxcheck[1:k,x] <- NA
    #  auxcheck <- na.omit(auxcheck)
    #  for(j in 1:length(x)) {
    #    if(length(unique(auxcheck[,x[j]]))<2) {
    #      x2del <- c(x2del,x[j])
    #      end <- 1
    #      }
    #    }
    #  if(end==1) break()
    #  }
    }
  x2del
  }

# perform ols estimation (internal use only)
doLS <- function(formula,group,xqual,data) {
  formOK <- formula
  if(length(xqual)>0) {
    auxform <- scanForm(formula)    
    xall <- names(auxform$lpar)
    form0 <- as.character(formula)
    xstr0 <- strsplit(gsub(" ","",form0[3]),"\\+")[[1]]
    nm0 <- max(extrLPar(formula)[,2])
    x2del <- checkQV(x=xqual,z=setdiff(xall,xqual),group=group,k=nm0+1,data=data)
    xstr0 <- setdiff(xstr0,x2del)
    formOK <- formula(paste(form0[2],"~",paste(xstr0,collapse="+"),sep=""))
    }
  res <- lm(formOK,data=data)
  res$call$formula <- formOK
  res
  }

# fit a distributed-lag linear regression model (internal use only)
dlaglm <- function(formula,group,xqual,data,adapt,no.select,max.gestation,min.width,max.lead,sign,selection,L) {
  auxscan <- scanForm(formula)
  y <- auxscan$y
  if(length(auxscan$X)>0) {
    lagPar <- auxscan$lpar
    lagType <- auxscan$ltype
    lagNam <- names(lagPar)  
    if(length(lagNam)>0) {
      if(is.null(group)) {
        maxlag <- nrow(na.omit(data[,c(y,lagNam)]))-2
        } else {
        auxmlag <- c()
        gruppi <- levels(factor(data[,group]))
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
      formOK <- creatForm(y,names(lagPar),group,lagType,lagPar)
      modOK <- doLS(formula=formOK,group=group,xqual=xqual,data=data)
      } else {
      bestPar <- vector("list",length=length(lagNam))
      names(bestPar) <- lagNam
      xOK <- no.select  #####
      fine <- 0
      while(fine==0) {
        xtest <- setdiff(lagNam,xOK)
        if(length(xtest)>0) {
          currentAIC <- rep(NA,length(xtest)) 
          names(currentAIC) <- xtest
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
              #####
              aux0 <- which(gammaCons$sx==0,arr.ind=T)
              aux1 <- which(gammaCons$dx[aux0]<=maxld,arr.ind=T)
              aux2 <- which.max(gammaCons$dx[aux0][aux1])
              gamPar0 <- as.numeric(c(rownames(gammaCons$dx)[aux0[aux2,1]],colnames(gammaCons$dx)[aux0[aux2,2]]))
              for(j in 1:length(lagPar)) {
                if(sum(is.na(lagPar[[j]]))==0) {
                  if(lagType[names(lagPar)[j]]=="gamma.lag") {
                    lagPar[[j]] <- gamPar0
                    } else {
                    lagPar[[j]] <- c(0,maxld)  
                    }
                  }
                }
              #####
              auxcons <- searchGrid(maxgs,minwd,maxld,lagType[xtest[i]])
              aic0 <- bhat0 <- c()
              for(j in 1:nrow(auxcons)) {
                testType <- lagType[c(xOK,xtest[i])]     
                testPar <- lagPar[c(xOK,xtest[i])]                      
                testPar[[xtest[i]]] <- auxcons[j,]
                form0 <- creatForm(y,names(testPar),group,testType,testPar)
                mod0 <- doLS(formula=form0,group=group,xqual=xqual,data=data)
                summ0 <- summary(mod0)$coefficients
                ixall <- rownames(summ0)
                iauxlab <- sapply(ixall,extrName)
                ixlab <- ixall[which(iauxlab==xtest[i])]
                if(ixlab %in% rownames(summ0)) {     
                  bhat0[j] <- summ0[ixlab,1]
                  } else {
                  bhat0[j] <- NA
                  }
                if(selection=="aic") {
                  aic0[j] <- AIC(mod0)
                  } else {
                  aic0[j] <- BIC(mod0)
                  }
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
              currentAIC[xtest[i]] <- aic0[auxbest]
              } else {
              testType <- lagType[c(xOK,xtest[i])]     
              testPar <- lagPar[c(xOK,xtest[i])]
              form0 <- creatForm(y,names(testPar),group,testType,testPar)
              mod0 <- doLS(formula=form0,group=group,xqual=xqual,data=data)
              auxsumm <- summary(mod0)$coefficients
              if(xtest[i] %in% rownames(auxsumm)) {
                if(selection=="aic") {
                  currentAIC[xtest[i]] <- AIC(mod0)
                  } else {
                  currentAIC[xtest[i]] <- BIC(mod0)  
                  }
                } else {
                currentAIC[xtest[i]] <- Inf     
                }
              }
            }
          xOK <- c(xOK,names(currentAIC)[which.min(currentAIC)])
          } else {
          fine <- 1
          }
        }
      formOK <- creatForm(y,names(bestPar),group,lagType,bestPar)
      modOK <- doLS(formula=formOK,group=group,xqual=xqual,data=data)
      }
    } else {
    if(is.null(group)) {
      formOK <- formula(paste(y,"~1",sep=""))
      } else {
      formOK <- formula(paste(y,"~-1+",group,sep=""))        
      }
    modOK <- doLS(formula=formOK,group=group,xqual=xqual,data=data)
    }
  #if(L>0) {
  #  res0 <- residuals(modOK)
  #  W <- matrix(nrow=length(res0),ncol=length(res0))
  #  for(i in 1:length(res0)) {
  #    for(j in 1:length(res0)) {
  #      ijlag <- abs(xtime[as.numeric(names(res0)[i])]-xtime[as.numeric(names(res0)[j])])
  #      W[i,j] <- max(0,(1-ijlag/(1+L)))*res0[i]*res0[j]
  #      }
  #    }
  #  Xmat <- model.matrix(modOK)
  #  Imat <- solve(t(Xmat)%*%Xmat)
  #  modOK$vcov <- Imat%*%t(Xmat)%*%W%*%Xmat%*%Imat
  #  } else {
  #  modOK$vcov <- bvar*summary(modOK)$sigma^2
  #  }
  modOK
  }

# compute lag effects of a covariate (internal use only)
lagEff <- function(model,x,cumul,conf,lag) {
  formstr <- strsplit(gsub(" ","",as.character(model$call$formula)[3]),"\\+")[[1]]
  xall <- names(model$coefficients)
  auxlab <- sapply(xall,extrName)
  xlab <- xall[which(auxlab==x)]
  auxscan <- scanForm(model$call)
  if(auxscan$ltype[x]=="quec.lag") {
    sx <- auxscan$lpar[[x]][1]
    dx <- auxscan$lpar[[x]][2]
    imu <- model$coefficients[xlab]
    icov <- vcov(model)[xlab,xlab]
    if(is.null(lag)) {
      xgrid <- 0:dx
      } else {
      xgrid <- lag
      }
    iH <- matrix(lagwei(c(sx,dx),xgrid,"quec.lag"),ncol=1)
    } else if(auxscan$ltype[x]=="qdec.lag") {
    sx <- auxscan$lpar[[x]][1]
    dx <- auxscan$lpar[[x]][2]
    imu <- model$coefficients[xlab]
    icov <- vcov(model)[xlab,xlab]
    if(is.null(lag)) {
      xgrid <- 0:dx
      } else {
      xgrid <- lag
      }
    iH <- matrix(lagwei(c(sx,dx),xgrid,"qdec.lag"),ncol=1)
    } else if(auxscan$ltype[x]=="gamma.lag") {
    delta <- auxscan$lpar[[x]][1]
    lambda <- auxscan$lpar[[x]][2]
    imu <- model$coefficients[xlab]
    icov <- vcov(model)[xlab,xlab]
    if(is.null(lag)) {
      xa <- 0  #####
      maxgam <- qgamma(0.99,1/(1-delta),-log(lambda))+xa
      xgrid <- 0:maxgam
      } else {
      xgrid <- lag
      }
    iH <- matrix(lagwei(c(delta,lambda),xgrid,"gamma.lag"),ncol=1)
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
applyDiff <- function(x,group,time=NULL,data,k) {
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
    gruppi <- levels(factor(data[,group]))
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
  if(!is.null(group) && length(group)!=1) stop("Argument 'group' must contain a single variable name",call.=F)
  if(!is.null(group) && is.na(group)) group <- NULL
  if(!is.null(group)) {
    data[,group] <- factor(data[,group])
    gruppi <- levels(data[,group])
    if(length(gruppi)<2) stop("The group factor must have at least 2 unique values",call.=F)
    g.id <- as.numeric(data[,group])
    glab <- gruppi
    } else {
    g.id <- glab <- rep(1,nrow(data))  
    }
  auxvar <- setdiff(c(x,group,time),colnames(data))
  if(length(auxvar)>0) stop("Unknown variable: '",auxvar[1],"'",sep="",call.=F)
  for(i in 1:length(x)) {
    if(isQuant(data[,x[i]])==F) stop("'",x[i],"' is not a numerical variable",call.=F)
    }
  if(!is.null(time) && is.na(time)) time <- NULL
  if(!is.null(time) && length(time)!=1) stop("Argument 'time' must contain a single variable name",call.=F)
  if(!is.null(time) && isQuant(data[,time])==F) stop("'time' must be a numerical variable",call.=F)
  if(length(combine)!=1 || (combine %in% c("choi","demetrescu"))==F) stop("Argument 'combine' must be either 'choi' or 'demetrescu'",call.=F)
  if(length(k)!=1 || !is.numeric(k) || k<0 || k!=round(k)) stop("Argument 'k' must be an non-negative integer",call.=F)
  if(length(lshort)!=1 || !is.logical(lshort)) stop("Argument 'lshort' must be a logical value",call.=F)  
  data[which(abs(data[,x])==Inf),x] <- NA
  if(!is.null(time)) {
    for(i in 1:length(g.id)) {
      auxind <- which(g.id==gruppi[i])
      idat <- data[auxind,]
      data[auxind,] <- idat[order(idat[,time]),]
      }
    }
  res <- list()
  for(i in 1:length(x)) {
    res[[i]] <- doADFtest(data[,x[i]],g.id=g.id,test=test,combine=combine,k=k,lshort=lshort,glab=glab)
    }
  names(res) <- x
  res
  }

# unit root test for a variable (internal use only)  
doADFtest <- function(x,g.id,test,combine,k,lshort,glab) {
  gruppi <- sort(unique(g.id))
  ifelse(test=="adf", auxalt <- "stationary", auxalt <- "unit root")
  auxstat <- auxp <- nwm <- c()
  options(warn=-1)
  for(i in 1:length(gruppi)) {
    auxind <- which(g.id==gruppi[i])
    auxdat <- na.omit(x[auxind])
    nwm[i] <- length(auxdat) #sum(!is.na(auxdat))
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
      } else {
      auxstat[i] <- NA
      auxp[i] <- NA
      }
    }
  if(length(auxstat)>1) names(auxstat) <- names(nwm) <- glab
  if(length(auxp)>0) {
    auxpStar <- auxp[which(auxp<1)]
    if(length(auxpStar)>0) {
      m <- length(auxpStar)
      logp <- qnorm(auxpStar)
      rhat <- 1-var(logp)
      rstar <- max(rhat,-1/(m-1))
      if(combine=="demetrescu") {
        auxz <- sum(logp)/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
        } else if(combine=="choi") {
        auxz <- sum(logp)/sqrt(m)
        }
      auxpval <- 2*pnorm(-abs(auxz))
      } else {
      auxz <- NA
      auxpval <- NA
      }
    res <- list(statistic=auxstat,n=nwm,z.value=auxz,p.value=auxpval)
    } else {
    res <- list(statistic=NULL,n=nwm,z.value=NULL,p.value=NULL)
    }
  options(warn=0)
  res
  }

# estimate parameters of the imputation model (internal use only)
impuFit <- function(xcont,xqual,group,data) {
  res <- G <- list()
  logL <- 0
  z.names <- c(group,xqual)
  for(i in 1:length(xcont)) {
    if(i==1) {
      if(is.null(z.names)) ipar <- "1" else ipar <- z.names  
      } else {
      ipar <- c(xcont[1:(i-1)],z.names)        
      }
    iform <- formula(paste(xcont[i],"~",paste(ipar,collapse="+"),sep=""))
    imod <- lm(iform,data=data)
    imod$call$formula <- iform
    #imod <- doLS(formula=iform,group=group,xqual=xqual,data=data)
    logL <- logL-0.5*AIC(imod,k=0)
    res[[i]] <- imod
    G[[i]] <- ipar
    }
  names(res) <- names(G) <- xcont
  list(hat=res,G=G,logL=logL)
  }

# predict missing values from the imputation model (internal use only)
impuPred <- function(mod,data) {
  est <- mod$hat
  G <- mod$G
  res <- data
  options(warn=-1)
  for(i in 1:length(est)) {
    iest <- est[[i]]
    inam <- names(est)[i]
    ipar <- G[[inam]]
    ina <- which(is.na(data[,inam]))
    if(length(ina)>0) res[ina,inam] <- predict(iest,res[ina,ipar])
    }
  options(warn=0)
  res
  }

# linear interpolation (internal use only)
linImp <- function(x) {
  res <- x
  auxNA <- which(is.na(x))
  if(length(auxNA)>0) {
    naL <- split(auxNA,cumsum(c(1,diff(auxNA)!=1)))
    for(i in 1:length(naL)) {
      ina <- naL[[i]]
      x1 <- min(ina)-1
      x2 <- max(ina)+1
      y1 <- x[x1]
      y2 <- x[x2]
      b <- (y2-y1)/(x2-x1)
      a <- y1-b*x1
      res[ina] <- a+b*ina
      }
    }
  res
  }

# spline reconstruction (internal use only)
splinImp <- function(x,group,data) {
  #
  doImp <- function(x) {
    res <- x
    auxNA <- which(is.na(x))
    if(length(auxNA)>0) {
      auxO <- which(!is.na(x))
      auxI <- intersect(auxNA,min(auxO):max(auxO))
      yI <- spline(1:length(x),x,xout=1:length(x))
      res[auxI] <- yI$y[auxI]
      }
    res
    }
  #
  if(!is.null(group)) {
    data[,group] <- factor(data[,group])
    gruppi <- levels(data[,group])    
    res <- data
    for(w in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[w])
      for(i in 1:length(x)) {
        res[auxind,x[i]] <- doImp(data[auxind,x[i]])
        }
      }    
    } else {
    res <- data
    for(i in 1:length(x)) {
      res[,x[i]] <- doImp(data[,x[i]])
      }
    }
  res
  }

# imputation of missing data (internal use only)
EM.imputation <- function(xcont,xqual,group,data,tol,maxiter,quiet=F) {
  #isNA <- matrix(0,nrow(data),ncol=length(xcont))
  #colnames(isNA) <- xcont
  #currentDat <- data
  #for(i in 1:length(xcont)) {
  #  ina <- which(is.na(data[,xcont[i]]))
  #   if(length(ina)>0) {
  #    isNA[ina,xcont[i]] <- 1
  #    #currentDat[ina,xcont[i]] <- 0
  #    }
  #  }
  nmiss <- apply(data[,xcont],2,function(v){sum(is.na(v))})
  xcont <- xcont[order(nmiss)]
  currentDat <- data
  for(i in 1:length(xcont)) {
    currentDat[which(is.na(currentDat[,xcont[i]])),xcont[i]] <- mean(data[,xcont[i]],na.rm=T)
    }
  currentFit <- impuFit(xcont=xcont,xqual=xqual,group=group,data=currentDat)
  currentLik <- -Inf
  #currentLik <- currentFit$logL
  fine <- forcend <- 0
  count <- 1
  if(quiet==F) {
    cat("Starting EM...")
    flush.console() 
    }
  while(fine==0) {
    newDat <- impuPred(currentFit,data)
    newFit <- impuFit(xcont=xcont,xqual=xqual,group=group,data=newDat)
    newLik <- newFit$logL
    if(quiet==F) {
      cat('\r',"EM iteration ",count,". Log-likelihood: ",newLik,sep="")
      flush.console() 
      }
    if(newLik<currentLik) {
      newLik <- currentLik
      newDat <- currentDat
      #warning("Forced stop of EM algorithm because likelihood has decreased",call.=F)
      fine <- 1
      } else {
      if(newLik-currentLik>tol & count<maxiter) {
        currentFit <- newFit
        currentLik <- newLik
        currentDat <- newDat
        count <- count+1
        } else {
        fine <- 1
        if(count>=maxiter) forcend <- 1
        }
      }
    }
  if(quiet==F) {
    if(forcend==0) {
      cat('\r',"EM converged after ",count," iterations. Log-likelihood: ",newLik,sep="","\n")
      } else {
      cat('\r',"EM stopped after ",maxiter," iterations. Log-likelihood: ",newLik,sep="","\n")      
      }
    }
  newDat
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
  auxs <- which(bmat[,1]!=0)
  bmat_s <- bmat[c(max(1,min(auxs)-1):min(nrow(bmat),max(auxs)+1)),]
  bval <- as.numeric(rownames(bmat_s))
  xgrid <- sort(unique(c(bval,seq(min(bval),max(bval),length=100))))
  ygrid <- cbind(spline(bval,bmat_s[,1],xout=xgrid)$y,spline(bval,bmat_s[,2],xout=xgrid)$y,spline(bval,bmat_s[,3],xout=xgrid)$y)
  #####
  auxsgn <- sign(bmat[auxs[1],1])
  if(auxsgn==1) {
    auxdel <- which(ygrid[,1]<0)
    } else {
    auxdel <- which(ygrid[,1]>0)
    }
  if(length(auxdel)>0) {
    auxInt <- c()
    for(i in 1:length(bval)) {
      auxInt[i] <- which(xgrid==bval[i])
      }
    for(i in 1:length(auxdel)) {
      isx <- auxInt[max(which(auxInt<=auxdel[i]))]
      idx <- auxInt[min(which(auxInt>=auxdel[i]))]
      ygrid[isx:idx,] <- NA
      }
    ygrid[c(1,nrow(ygrid)),] <- 0
    for(j in 1:ncol(ygrid)) {
      ygrid[,j] <- linImp(ygrid[,j])
      }
    }
  #####
  plot(0,type="n",xlim=c(min(xaux),max(xaux)),ylim=c(lowLim,upLim),yaxs="i",xaxs="i",cex.lab=1.2,
    lwd=2,xaxt="n",yaxt="n",xlab="Lag",ylab="Coefficient",main=title,cex.main=1.2) 
  if(cumul==T) {
    mtext("cumulative lag shape",cex=0.9)
    } else {
    mtext("instantaneous lag shape",cex=0.9)    
    }
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

# check missing values (internal use only)
checkNA <- function(x,group,data) {
  if(sum(!is.na(data[,x]))<3) stop("Variable '",x,"' has less than 3 observed values",call.=F)
  if(!is.null(group)) {
    gruppi <- levels(factor(data[,group]))
    for(i in 1:length(gruppi)) {
      auxind <- which(data[,group]==gruppi[i])
      if(sum(!is.na(data[auxind,x]))<1) {
        stop("Variable '",x,"' has no observed values in group '",gruppi[i],"'",call.=F)
        }
      }
    }
  }

# fit a dlsem
dlsem <- function(model.code,group=NULL,time=NULL,exogenous=NULL,data,log=FALSE,diff.options=list(test="adf",combine="choi",k=0,lshort=TRUE,maxdiff=3),
  imput.options=list(tol=0.0001,maxiter=500,recons=TRUE,no.imput=NULL),
  global.control=NULL,local.control=NULL) {
  #
  if(!is.list(model.code) || sum(sapply(model.code,class)!="formula")>0) stop("Argument 'model code' must be a list of formulas",call.=F)
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be an object of class 'data.frame'",call.=F)  
  if(nrow(data)<3) stop("Data must contain at least 3 records",call.=F)
  nameOfData <- deparse(substitute(data))
  if(!is.null(group) && length(group)!=1) stop("Argument 'group' must contain a single variable name",call.=F)
  if(!is.null(group) && is.na(group)) group <- NULL
  if(!is.null(group)) {
    if(length(group)!=1) stop("Argument 'group' must be of length 1",call.=F)
    if((group %in% colnames(data))==F) stop("Unknown variable '",group,"' provided to argument 'group'",call.=F)
    gruppi <- levels(factor(data[,group]))
    if(length(gruppi)<2) stop("The group factor must have at least 2 unique values",call.=F)
    data[,group] <- factor(data[,group])
    }
  if(!is.null(time)) {
    if(length(time)!=1) stop("Argument 'time' must be of length 1",call.=F)
    if((time %in% colnames(data))==F) stop("Unknown variable '",time,"' provided to argument 'time'",call.=F)
    if(!is.null(group)) {
      for(i in 1:length(gruppi)) {
        iind <- which(data[,group]==gruppi[i])
        idat <- data[iind,]
        data[iind,] <- idat[order(idat[,time]),]
        }
      } else {
      data <- data[order(data[,time]),]
      }
    }
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
  if(is.null(imput.options)) imput.options <- list(tol=0.0001,maxiter=500,recons=T,no.imput=NULL)
  auxch <- setdiff(names(imput.options),c("tol","maxiter","recons","no.imput"))
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
  if("recons" %in% names(imput.options)) {
    recons <- imput.options$recons
    } else {
    recons <- T
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
  if((global.control$selection %in% c("aic","bic"))==F) stop("Component 'selection' in argument 'global.options' must be one among 'aic' and 'bic'",call.=F)
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
        if(icheck$ltype[j]=="quec.lag") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>icheck$lpar[[j]][2]) {
            stop("Invalid settings for lag shape 'quec.lag' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="qdec.lag") {
          if(length(icheck$lpar[[j]])!=2 || !identical(icheck$lpar[[j]],round(icheck$lpar[[j]])) || sum(icheck$lpar[[j]]<0)>0 || icheck$lpar[[j]][1]>icheck$lpar[[j]][2]) {
            stop("Invalid settings for lag shape 'qdec.lag' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
            }
          } else if(icheck$ltype[j]=="gamma.lag") {
          if(length(icheck$lpar[[j]])!=2 || sum(icheck$lpar[[j]]<=0)>0 || sum(icheck$lpar[[j]]>=1)>0) {
            stop("Invalid settings for lag shape 'gamma.lag' in '",paste(icheck$y," ~ ",paste(icheck$X,collapse="+"),sep=""),"'",call.=F)
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
    }
  for(i in 1:length(nodenam)) {
    if(isQuant(data[,nodenam[i]])==F) stop("Qualitative variables cannot appear in 'model.code': ",nodenam[i],sep="",call.=F)
    }
  if(!is.null(exogenous)) {
    for(i in 1:length(exogenous)) {
      if(isQuant(data[,exogenous[i]])==F) {
        if(sum(is.na(data[,exogenous[i]]))>0) stop("Qualitative variables cannot contain missing values: ",exogenous[i],sep="",call.=F)
        }
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
  if(is.null(topG)) stop("The DAG contains directed cycles",call.=F)  
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
  xfact <- c()
  for(i in 1:length(nodenam)) {
    if(isQuant(data[,nodenam[i]])==F) {
      xfact <- c(xfact,nodenam[i])
      data[,nodenam[i]] <- factor(data[,nodenam[i]])
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
  if(recons==T) {
    x2recons <- setdiff(difftest,no.imput)
    if(length(x2recons)>0) data <- splinImp(x2recons,group,data)
    }
  ndiff <- 0
  if(length(difftest)>0 & maxdiff>0) {
    fine <- 0
    cat("Checking stationarity...")
    flush.console()
    while(fine==0) {
      auxp <- c()
      urtList <- unirootTest(difftest,group=group,data=data,test=test,combine=combine,k=k,lshort=lshort)
      for(i in 1:length(difftest)) {
        ipvl <- urtList[[i]]$p.value
        if(is.null(ipvl)) {
          auxp[i] <- 0
          } else {
          if(is.na(ipvl)) {
            auxp[i] <- 0
            } else {
            auxp[i] <- ipvl
            }
          }
        }
      if(sum(auxp>0.05)>0) {
        if(ndiff<maxdiff) {
          ndiff <- ndiff+1
          data <- applyDiff(x=difftest,group=group,data=data,k=rep(ndiff,length(difftest)))                                        
          } else {
          fine <- 1
          }
        } else {
        fine <- 1
        }
      }
    } else {
    urtList <- NULL  
    }
  if(ndiff==0) {
    cat('\r',"No differentiation performed")
    } else {
    cat('\r',"Order",ndiff,"differentiation performed")    
    }
  cat("\n")
  if(is.null(group)) {
    nK <- length(nodenam)+1
    } else {
    nK <- length(nodenam)+nlevels(factor(data[,group]))
    }
  auxna <- apply(data[,nodenam],1,function(x){sum(is.na(x))})
  if(sum(auxna)>0) {
    auxOK <- unname(which(auxna<length(nodenam))) 
    if(ndiff>0) {
      if(is.null(group)) {
        auxind <- 1:ndiff
        } else {
        auxind <- c()
        for(i in 1:length(gruppi)) {
          auxind <- c(auxind,which(data[,group]==gruppi[i])[1:ndiff])          
          }
        }                              
      auxOK <- setdiff(auxOK,auxind)                    
      }
    if(sum(is.na(data[auxOK,]))>0 && maxiter>0) {
      x2imp <- setdiff(nodenam,c(no.imput,xfact))
      if(length(x2imp)>0) {
        for(i in 1:length(x2imp)) {
          nachk <- checkNA(x2imp[i],group,data[auxOK,])
          }
        nIm <- sum(is.na(data[,x2imp]))
        if(nIm>0) {
          data[auxOK,c(group,xfact,x2imp)] <- EM.imputation(xcont=x2imp,xqual=xfact,group=group,data=data[auxOK,c(group,xfact,x2imp),drop=F],tol=tol,maxiter=maxiter)
          }
        }
      }
    }
  nomi <- c()
  optList <- list()
  cat("Starting estimation...")
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
      iform <- formula(paste(as.character(model.code[[i]])[2],"~",paste(exogenous,collapse="+"),"+",as.character(model.code[[i]])[3],sep=""))
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
    res[[i]] <- dlaglm(iform,group=group,xqual=xfact,data=data,adapt=iad,no.select=exogenous,max.gestation=iges,min.width=iwd,max.lead=ilead,sign=isg,selection=global.control$selection)  ### L=iL,
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
  callList <- lapply(res,function(z){z$call})
  out <- list(estimate=res[topG],call=callList,exogenous=exogenous,group=group,log=log,uniroot=urtList,ndiff=ndiff,
    diff.options=diff.options,imput.options=imput.options,selection=global.control$selection,adaptation=optList,
    Rsq=RsqCalc(res[topG]),data.orig=origdat[,c(group,time,nodenam)],data.used=data[,c(group,time,nodenam)])
  class(out) <- "dlsem"
  out
  }

# automated plots of lag shapes
auto.lagPlot <- function(x,cumul=FALSE,conf=0.95,plotDir=NULL) {
  if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  if(length(cumul)!=1 || !is.logical(cumul)) stop("Arguent 'cumul' must be a logical value",call.=F)
  if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  if(is.null(plotDir)) plotDir <- getwd()
  for(i in 1:length(x$estimate)) {  
    scan0 <- scanForm(x$estimate[[i]]$call$formula)
    inam <- scan0$y
    isumm <- summary(x$estimate[[inam]])$coefficients
    ilagged <- names(scan0$ltype)[which(scan0$ltype!="none")]
    ilab <- sapply(rownames(isumm),extrName)
    if(length(ilagged)>0) {                                  
      for(j in 1:length(ilagged)) {
        ijlab <- names(ilab)[which(ilab==ilagged[j])]
        if(isumm[ijlab,4]<=1-conf) {
          pdf(file.path(plotDir,paste(inam,"~",ilagged[j],".pdf",sep="")))
          lagPlot(x,path=paste(ilagged[j],"*",inam,sep=""),cumul=cumul,conf=conf)
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
    cat(" Number of groups: ",nlevels(x$data.used[,x$group]),"\n",sep="")
    } else {
    cat(" No groups","\n")
    }
  cat(" Number of endogenous variables: ",length(x$estimate),"\n",sep="")
  if(!is.null(x$exogenous)) {
    cat(" Number of exogenous variables: ",length(x$exogenous),"\n",sep="")
    } else {
    cat(" No exogenous variables","\n")
    }
  if(n.e>0) {
    cat(" ",n.e,"/",N.e," significant edges at 5% level","\n",sep="")
    } else {
    cat(" No significant edges at 5% level","\n")  
    }
  }

# format summary table (internal use only)
formatSumm <- function(x,newnam) {
  if(newnam==T) {
    colnam <- c("theta","se(theta)","t value","Pr(>|t|)","")
    } else {
    colnam <- c(colnames(x),"")
    }
  auxp <- rep("",nrow(x))
  for(i in 1:nrow(x)) {
    ip <- x[i,ncol(x)]
    if(ip<0.1) auxp[i] <- "."
    if(ip<0.05) auxp[i] <- "*"
    if(ip<0.01) auxp[i] <- "**"
    if(ip<0.001) auxp[i] <- "***"
    }
  res <- data.frame(x,auxp)
  colnames(res) <- colnam
  res
  }

# summary method for class dlsem
summary.dlsem <- function(object,...) {
  elev <- glev <- c()
  enam <- object$exogenous 
  if(!is.null(enam)) {
    for(i in 1:length(enam)) {
      if(!is.factor(object$data.used[,enam[i]])) {
        elev <- c(elev,enam[i]) 
        } else {
        elev <- c(elev,paste(enam[i],levels(object$data.used[,enam[i]]),sep=""))
        }
      }
    }
  if(!is.null(object$group)) {
    glev <- paste(object$group,levels(object$data.used[,object$group]),sep="")
    }
  estim <- object$estimate
  fitI <- c()
  fitI[1] <- RsqCalc(estim)["(overall)"]
  fitI[2] <- AIC(object)["(overall)"]
  fitI[3] <- BIC(object)["(overall)"]
  names(fitI) <- c("Rsq","AIC","BIC")
  summList <- summList_e <- summList_g <- vector("list",length=length(estim))
  names(summList) <- names(summList_e) <- names(summList_g) <- names(estim)
  summS <- matrix(nrow=length(estim),ncol=2)
  rownames(summS) <- names(estim)
  colnames(summS) <- c("sigma","df")
  for(i in 1:length(estim)) {
    iB <- summary(estim[[i]])$coefficients
    inam <- setdiff(rownames(iB),c("(Intercept)",elev,glev))
    if(length(inam)>0) {
      summList[[i]] <-formatSumm(iB[inam,,drop=F],newnam=F)
      }
    if(length(elev)>0) {
      elevOK <- intersect(rownames(iB),elev)
      if(length(elevOK)>0) {
        summList_e[[i]] <- formatSumm(iB[elevOK,,drop=F],newnam=F)
        } else {
        summList_e[[i]] <- NULL
        }
      }
    if(length(glev)>0) {
      summList_g[[i]] <- formatSumm(iB[glev,,drop=F],newnam=F)   
      } else {
      summList_g[[i]] <- formatSumm(iB["(Intercept)",,drop=F],newnam=F)          
      }
    summS[i,] <- c(summary(estim[[i]])$sigma,df=estim[[i]]$df.residual)
    }
  OUT <- list(endogenous=summList,exogenous=summList_e,group=summList_g,sigma=summS,gof=fitI)
  class(OUT) <- "summary.dlsem"
  OUT
  }

# print method for class summary.dlsem
print.summary.dlsem <- function(x,...) {
  cat("ENDOGENOUS PART","\n","\n")
  for(i in 1:length(x$endogenous)) {
    cat("Response: ",names(x$endogenous)[i],sep="","\n")
    isumm <- x$endogenous[[i]]
    if(!is.null(isumm)) {
      print(isumm)
      } else {
      cat("-","\n")
      }
    if(i<length(x$endogenous)) cat("\n")
    }
  cat("\n","\n")
  cat("EXOGENOUS PART","\n")
  if(sum(sapply(x$exogenous,is.null))==0) {
    cat("\n")
    for(i in 1:length(x$exogenous)) {
      cat("Response: ",names(x$exogenous)[i],sep="","\n")
      print(x$exogenous[[i]])
      cat("\n")
      }
    } else {
    cat(" -","\n","\n")
    }
  fitI <- x$gof
  cat("\n")
  cat("INTERCEPTS","\n")
  if(sum(sapply(x$group,is.null))==0) {
    cat("\n")
    for(i in 1:length(x$group)) {
      cat("Response: ",names(x$group)[i],sep="","\n")
      print(x$group[[i]])
      cat("\n")
      }
    } else {
    cat(" -","\n","\n")
    }
  fitI <- x$gof
  cat("\n")
  cat("GOODNESS OF FIT","\n","\n")
  cat("R-squared: ",round(fitI[1],4),"\n",sep="")
  cat("AIC: ",fitI[2],"\n",sep="")
  cat("BIC: ",fitI[3],"\n",sep="")
  }

# nobs method for class dlsem
nobs.dlsem <- function(object,...) {
  lapply(object$estimate,nobs)
  }

# coef method for class dlsem
coef.dlsem <- function(object,...) {
  lapply(object$estimate,coef)
  }

# vcov method for class dlsem
vcov.dlsem <- function(object,...) {
  lapply(object$estimate,vcov)
  }

# logLik method for class dlsem
logLik.dlsem <- function(object,...) {
  lapply(object$estimate,logLik)
  }

# AIC method for class dlsem
AIC.dlsem <- function(object,...) {
  res <- sapply(object$estimate,AIC)
  OUT <- c(res,sum(res))
  names(OUT) <- c(names(res),"(overall)")
  OUT
  }

# BIC method for class dlsem
BIC.dlsem <- function(object,...) {
  res <- sapply(object$estimate,BIC)
  OUT <- c(res,sum(res))
  names(OUT) <- c(names(res),"(overall)")
  OUT
  }

# format fitted or residuals of an object of class dlsem (internal use only)
formatFit <- function(object,pred) {
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

# fitted method for class dlsem
fitted.dlsem <- function(object,...) {
  pred <- list()
  for(i in 1:length(object$estimate)) {
    pred[[i]] <- fitted(object$estimate[[i]],...)
    }
  names(pred) <- names(object$estimate)
  formatFit(object,pred)
  }

# residuals method for class dlsem
residuals.dlsem <- function(object,...) {
  pred <- list()
  for(i in 1:length(object$estimate)) {
    pred[[i]] <- residuals(object$estimate[[i]],...)
    }
  names(pred) <- names(object$estimate)
  formatFit(object,pred)
  }

# predict method for class dlsem
predict.dlsem <- function(object,newdata=NULL,...) {  
  pred <- list()
  if(is.null(newdata)) {
    for(i in 1:length(object$estimate)) {
      pred[[i]] <- predict(object$estimate[[i]],...)
      }
    if(!is.null(object$group)) {
      igrou <- as.character(object$data.used[,object$group])
      names(igrou) <- rownames(object$data.used)
      pred <- c(list(igrou),pred)
      }
    } else {
    Z <- newdata
    if(ncol(Z)>0) {
      Xq <- c()
      for(i in 1:ncol(Z)) {
        if(isQuant(Z[,i])) Xq <- c(Xq,colnames(Z)[i])
        }
      if(length(Xq)>0) {
        if(object$log==T) {
          for(i in 1:length(Xq)) {
            if(sum(Z[,Xq[i]]<=0,na.rm=T)==0) Z[,Xq[i]] <- log(Z[,Xq[i]])
            }
          }
        if(object$ndiff>0) Z <- applyDiff(x=Xq,group=object$group,data=Z,k=rep(object$ndiff,length(Xq))) 
        }
      }
    pred <- list()
    for(i in 1:length(object$estimate)) {
      pred[[i]] <- predict(object$estimate[[i]],newdata=Z,...)
      }
    if(!is.null(object$group)) {
      igrou <- as.character(Z[,object$group])  ####
      names(igrou) <- rownames(Z)
      pred <- c(list(igrou),pred)
      }    
    }
  names(pred) <- c(object$group,names(object$estimate))
  auxind <- lapply(pred,names)
  ind <- sort(unique(as.numeric(unlist(auxind))))
  res <- data.frame(matrix(nrow=length(ind),ncol=length(pred)))
  colnames(res) <- names(pred)
  rownames(res) <- ind
  for(i in 1:length(pred)) {
    inam <- names(pred[[i]])
    for(j in 1:length(inam)) {
      res[inam[j],i] <- pred[[i]][inam[j]]
      }
    }
  if(!is.null(object$group)) {
    res[,object$group] <- factor(res[,object$group],levels=levels(object$data.used[,object$group]))
    }
  res
  }

# compute R-squared (internal use only)
RsqCalc <- function(xfit) {
  expdev <- resdev <- Rsq <- c()
  for(i in 1:length(xfit)) {
    auxss <- anova(xfit[[i]])$'Sum Sq'
    expdev[i] <- sum(auxss[1:(length(auxss)-1)])
    resdev[i] <- auxss[length(auxss)]
    Rsq[i] <- summary(xfit[[i]])$'r.squared'
    }
  OUT <- c(Rsq,sum(expdev)/(sum(expdev)+sum(resdev)))
  names(OUT) <- c(names(xfit),"(overall)")
  OUT
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
      ijnam <- extrName(auxnam[j])
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

# compute the coefficient associated to each edge at different time lags (internal use only)
edgeCoeff <- function(x,lag=NULL,conf=0.95) {
  #if(class(x)!="dlsem") stop("Argument 'x' must be an object of class 'dlsem'",call.=F)
  #if(length(conf)!=1 || !is.numeric(conf) || conf<=0 || conf>=1) stop("Arguent 'conf' must be a real number greater than 0 and less than 1",call.=F)
  #if(!is.null(lag)) {
  #  for(i in 1:length(lag)) {
  #    if(!is.numeric(lag[i]) || is.na(lag[i]) || lag[i]<0 || lag[i]!=round(lag[i])) stop("Argument 'lag' must contain non-negative integer numbers only",call.=F)
  #    }  
  #  }
  nomi <- names(x$estimate)
  laglen <- c()
  for(i in 1:length(nomi)) {
    isumm <- summary(x$estimate[[nomi[i]]])$coefficients
    auxnam <- rownames(isumm)
    for(j in 1:length(auxnam)) {
      ijnam <- extrName(auxnam[j])
      if(ijnam %in% nomi) {
        cumb <- lagEff(model=x$estimate[[nomi[i]]],x=ijnam,cumul=F,conf=conf,lag=NULL)
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
      ijnam <- extrName(auxnam[j])
      if(ijnam %in% nomi) {
        for(w in 1:length(lagOK)) {
          bList[[w]] <- rbind(bList[[w]],lagEff(model=x$estimate[[nomi[i]]],x=ijnam,cumul=F,conf=conf,lag=lagOK[w]))
          rownames(bList[[w]])[nrow(bList[[w]])] <- paste(ijnam,"~",nomi[i],sep="")
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
  #
  #
  par(xpd=T)
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
  if(!is.null(var1) && length(var1)!=1) stop("Argument 'var1' must be of length 1",call.=F)
  if(is.null(var2) || is.na(var2)) stop("Argument 'var2' is missing",call.=F)  
  if(!is.null(var2) && length(var2)!=1) stop("Argument 'var2' must be of length 1",call.=F)
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
  if(!is.character(to)) stop("Invalid argument 'to'",call.=F)
  if(length(to)!=1) stop("Argument 'to' must be of length 1",call.=F)
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
        auxeff <- lagEff(model=x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf,lag=NULL)
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
    quan <- -qnorm((1-conf)/2)
    #
    sd_calc <- function(muvet,sdvet) { sqrt(prod(muvet^2+sdvet^2)-prod(muvet^2)) }
    #
    sd_sum <- function(sdvet) {
      res <- c()
      for(i in 1:length(sdvet)) {
        res[i] <- sqrt(sum(sdvet[1:i]^2))
        }
      res
      }    
    #
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
      outList[[i]] <- matrix(nrow=length(lagOK),ncol=2)
      rownames(outList[[i]]) <- lagOK
      colnames(outList[[i]]) <- c("estimate","std. error")
      auxbetalag <- list()
      for(j in 2:length(newPathList[[i]])) {
        auxnam <- paste(newPathList[[i]][j],"~",newPathList[[i]][j-1],sep="")
        auxeff <- lagEff(model=x$estimate[[newPathList[[i]][j]]],x=newPathList[[i]][j-1],cumul=F,conf=conf,lag=lagOK)
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
          outList[[i]][j,] <- c(mupath,sdpath)
          } else {
          outList[[i]][j,] <- rep(0,2)
          }
        }
      }
    names(outList) <- sapply(newPathList,function(x){paste(x,collapse="*")})
    out <- matrix(nrow=length(lagOK),ncol=2)
    rownames(out) <- lagOK
    colnames(out) <- c("estimate","std. error")
    for(j in 1:length(lagOK)) {
      #auxover <- rep(0,3)
      auxover <- rep(0,2)
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
          outList[[i]][,1] <- cumsum(outList[[i]][,1])
          outList[[i]][,2] <- sd_sum(outList[[i]][,2])
          }
        }    
      }
    #
    for(i in 1:length(outList)) {
      imu <- outList[[i]][,1]
      isd <- outList[[i]][,2]
      outList[[i]] <- cbind(imu,imu-quan*isd,imu+quan*isd)
      colnames(outList[[i]]) <- c("estimate",paste(c("lower ","upper "),conf*100,"%",sep=""))
      }
    #
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
