### R code from vignette source 'dlsem_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_vignette.Rnw:295-296 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_vignette.Rnw:305-306 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_vignette.Rnw:338-339
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_vignette.Rnw:344-346
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_vignette.Rnw:387-392
###################################################
mycode <- list(
  Job ~ 1,
  Consum~quec(Job,0,15),
  Pollution~quec(Job,0,15)+quec(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_vignette.Rnw:437-444
###################################################
mycontrol <- list(
  adapt=c(Consum=T,Pollution=T),
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  min.width=list(Consum=c(Job=5),Pollution=c(Job=5,Consum=5)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 7: dlsem_vignette.Rnw:486-488
###################################################
mod0 <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,control=mycontrol,log=T)


###################################################
### code chunk number 8: dlsem_vignette.Rnw:508-509 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 9: dlsem_vignette.Rnw:531-532
###################################################
summary(mod0)


###################################################
### code chunk number 10: dlsem_vignette.Rnw:544-545
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 11: dlsem_vignette.Rnw:579-580
###################################################
causalEff(mod0,from="Job",to="Pollution",lag=seq(0,20,by=5),cumul=T)


###################################################
### code chunk number 12: dlsem_vignette.Rnw:604-605
###################################################
causalEff(mod0,from="Job",to="Pollution",cumul=T)$overall


###################################################
### code chunk number 13: dlsem_vignette.Rnw:616-618 (eval = FALSE)
###################################################
## lagPlot(mod0,path="Job*Pollution")
## lagPlot(mod0,path="Job*Consum*Pollution")


###################################################
### code chunk number 14: dlsem_vignette.Rnw:625-626 (eval = FALSE)
###################################################
## lagPlot(mod0,from="Job",to="Pollution")


