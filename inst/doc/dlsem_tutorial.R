### R code from vignette source 'dlsem_tutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_tutorial.Rnw:274-275
###################################################
require(dlsem)


###################################################
### code chunk number 2: dlsem_tutorial.Rnw:280-282
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 3: dlsem_tutorial.Rnw:333-338
###################################################
mycode <- list(
  Job ~ 1,
  Consum~quec(Job,0,15),
  Pollution~quec(Job,0,15)+quec(Consum,0,15)
  )


###################################################
### code chunk number 4: dlsem_tutorial.Rnw:391-398
###################################################
mycontrol <- list(
  adapt=c(Consum=T,Pollution=T),
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  min.width=list(Consum=c(Job=5),Pollution=c(Job=5,Consum=5)),
  max.width=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 5: dlsem_tutorial.Rnw:441-443
###################################################
mod0 <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,control=mycontrol,uniroot.check=T,log=T)


###################################################
### code chunk number 6: dlsem_tutorial.Rnw:466-467 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 7: dlsem_tutorial.Rnw:489-490
###################################################
summary(mod0)


###################################################
### code chunk number 8: dlsem_tutorial.Rnw:502-503
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 9: dlsem_tutorial.Rnw:537-538
###################################################
causalEff(mod0,from="Job",to="Pollution",lag=seq(0,20,by=5),cumul=T)


###################################################
### code chunk number 10: dlsem_tutorial.Rnw:562-563
###################################################
causalEff(mod0,from="Job",to="Pollution",cumul=T)


###################################################
### code chunk number 11: dlsem_tutorial.Rnw:574-576 (eval = FALSE)
###################################################
## lagPlot(mod0,path="Job*Pollution")
## lagPlot(mod0,path="Job*Consum*Pollution")


###################################################
### code chunk number 12: dlsem_tutorial.Rnw:583-584 (eval = FALSE)
###################################################
## lagPlot(mod0,from="Job",to="Pollution")


