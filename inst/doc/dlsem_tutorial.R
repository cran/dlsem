### R code from vignette source 'dlsem_tutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_tutorial.Rnw:396-397
###################################################
require(dlsem)


###################################################
### code chunk number 2: dlsem_tutorial.Rnw:402-404
###################################################
data(agres)
summary(agres)


###################################################
### code chunk number 3: dlsem_tutorial.Rnw:448-453
###################################################
mycode <- list(
  GVA~quec(NPATENT,0,15),
  PPI~quec(NPATENT,0,15)+quec(GVA,0,15),
  ENTR_INCOME~quec(NPATENT,0,15)+quec(GVA,0,15)
  )


###################################################
### code chunk number 4: dlsem_tutorial.Rnw:503-512
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
### code chunk number 5: dlsem_tutorial.Rnw:555-557
###################################################
mod0 <- dlsem(mycode,group="COUNTRY",context=c("GDP","FARM_SIZE"),
  data=agres,control=mycontrol,uniroot.check=T,imputation=T,log=T)


###################################################
### code chunk number 6: dlsem_tutorial.Rnw:564-565 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 7: dlsem_tutorial.Rnw:594-595
###################################################
summary(mod0)


###################################################
### code chunk number 8: dlsem_tutorial.Rnw:607-608
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 9: dlsem_tutorial.Rnw:633-634
###################################################
pathAnal(mod0,from="NPATENT",to="ENTR_INCOME",lag=c(5,10,15,20,25),cumul=T)


###################################################
### code chunk number 10: dlsem_tutorial.Rnw:638-639
###################################################
pathAnal(mod0,from="NPATENT",to="PPI",lag=c(5,10,15,20,25),cumul=T)


###################################################
### code chunk number 11: dlsem_tutorial.Rnw:663-664 (eval = FALSE)
###################################################
## lagPlot(mod0,from="NPATENT",to="ENTR_INCOME")


###################################################
### code chunk number 12: dlsem_tutorial.Rnw:668-669 (eval = FALSE)
###################################################
## lagPlot(mod0,from="NPATENT",to="PPI")


