### R code from vignette source 'dlsem_tutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_tutorial.Rnw:275-276
###################################################
require(dlsem)


###################################################
### code chunk number 2: dlsem_tutorial.Rnw:281-283
###################################################
data(agres)
summary(agres)


###################################################
### code chunk number 3: dlsem_tutorial.Rnw:332-337
###################################################
mycode <- list(
  GVA~quec(NPATENT,0,15),
  PPI~quec(NPATENT,0,15)+quec(GVA,0,15),
  ENTR_INCOME~quec(NPATENT,0,15)+quec(GVA,0,15)
  )


###################################################
### code chunk number 4: dlsem_tutorial.Rnw:388-397
###################################################
mycontrol <- list(
  adapt=c(GVA=T,PPI=T,ENTR_INCOME=T),
  max.gestation=list(GVA=c(NPATENT=3),PPI=c(NPATENT=3,GVA=3),
    ENTR_INCOME=c(NPATENT=3,GVA=3)),
  min.width=list(GVA=c(NPATENT=5),PPI=c(NPATENT=5,GVA=5),
    ENTR_INCOME=c(NPATENT=5,GVA=5)),
  sign=list(GVA=c(NPATENT="+"),PPI=c(NPATENT="-",GVA="-"),
    ENTR_INCOME=c(NPATENT="+",GVA="+"))
  )


###################################################
### code chunk number 5: dlsem_tutorial.Rnw:442-444
###################################################
mod0 <- dlsem(mycode,group="COUNTRY",exogenous=c("GDP","FARM_SIZE"),
  data=agres,control=mycontrol,uniroot.check=T,imputation=T,log=T)


###################################################
### code chunk number 6: dlsem_tutorial.Rnw:459-460 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 7: dlsem_tutorial.Rnw:492-493
###################################################
summary(mod0)


###################################################
### code chunk number 8: dlsem_tutorial.Rnw:505-506
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 9: dlsem_tutorial.Rnw:544-545
###################################################
pathAnal(mod0,from="NPATENT",to="ENTR_INCOME",lag=c(5,10,15,20,25),cumul=T)


###################################################
### code chunk number 10: dlsem_tutorial.Rnw:549-550
###################################################
pathAnal(mod0,from="NPATENT",to="PPI",lag=c(5,10,15,20,25),cumul=T)


###################################################
### code chunk number 11: dlsem_tutorial.Rnw:572-573 (eval = FALSE)
###################################################
## lagPlot(mod0,from="NPATENT",to="ENTR_INCOME")


###################################################
### code chunk number 12: dlsem_tutorial.Rnw:577-578 (eval = FALSE)
###################################################
## lagPlot(mod0,from="NPATENT",to="PPI")


