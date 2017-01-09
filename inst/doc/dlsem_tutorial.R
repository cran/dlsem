### R code from vignette source 'dlsem_tutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_tutorial.Rnw:278-279 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_tutorial.Rnw:288-289 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_tutorial.Rnw:308-309
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_tutorial.Rnw:314-316
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_tutorial.Rnw:371-376
###################################################
mycode <- list(
  Job ~ 1,
  Consum~quec(Job,0,15),
  Pollution~quec(Job,0,15)+quec(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_tutorial.Rnw:421-428
###################################################
mycontrol <- list(
  adapt=c(Consum=T,Pollution=T),
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  min.width=list(Consum=c(Job=5),Pollution=c(Job=5,Consum=5)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 7: dlsem_tutorial.Rnw:470-472
###################################################
mod0 <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,control=mycontrol,log=T)


###################################################
### code chunk number 8: dlsem_tutorial.Rnw:492-493 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 9: dlsem_tutorial.Rnw:515-516
###################################################
summary(mod0)


###################################################
### code chunk number 10: dlsem_tutorial.Rnw:528-529
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 11: dlsem_tutorial.Rnw:563-564
###################################################
causalEff(mod0,from="Job",to="Pollution",lag=seq(0,20,by=5),cumul=T)


###################################################
### code chunk number 12: dlsem_tutorial.Rnw:588-589
###################################################
causalEff(mod0,from="Job",to="Pollution",cumul=T)$overall


###################################################
### code chunk number 13: dlsem_tutorial.Rnw:600-602 (eval = FALSE)
###################################################
## lagPlot(mod0,path="Job*Pollution")
## lagPlot(mod0,path="Job*Consum*Pollution")


###################################################
### code chunk number 14: dlsem_tutorial.Rnw:609-610 (eval = FALSE)
###################################################
## lagPlot(mod0,from="Job",to="Pollution")


