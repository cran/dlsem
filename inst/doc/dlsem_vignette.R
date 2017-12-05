### R code from vignette source 'dlsem_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_vignette.Rnw:302-303 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_vignette.Rnw:310-311 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_vignette.Rnw:343-344
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_vignette.Rnw:349-351
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_vignette.Rnw:392-397
###################################################
mycode <- list(
  Job ~ 1,
  Consum~quec(Job,0,15),
  Pollution~quec(Job,0,15)+quec(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_vignette.Rnw:429-433 (eval = FALSE)
###################################################
## list(
##   min.width=list(Pollution=c(Consum=5)),
##   sign=list(Pollution=c(Job="+",Consum="+"))
##   )


###################################################
### code chunk number 7: dlsem_vignette.Rnw:448-450
###################################################
mycon_G <- list(adapt=T,max.gestation=3,max.lead=15,min.width=5,sign="+")
mycon_L <- list()


###################################################
### code chunk number 8: dlsem_vignette.Rnw:457-465
###################################################
mycon_G <- list()
mycon_L <- list(
  adapt=c(Consum=T,Pollution=T),
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  min.width=list(Consum=c(Job=5),Pollution=c(Job=5,Consum=5)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 9: dlsem_vignette.Rnw:470-476
###################################################
mycon_G <- list(adapt=T,min.width=5)
mycon_L <- list(
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 10: dlsem_vignette.Rnw:517-519
###################################################
mod0 <- dlsem(mycode,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=mycon_G,local.control=mycon_L,log=T)


###################################################
### code chunk number 11: dlsem_vignette.Rnw:529-530 (eval = FALSE)
###################################################
## plot(mod0)


###################################################
### code chunk number 12: dlsem_vignette.Rnw:552-553
###################################################
summary(mod0)


###################################################
### code chunk number 13: dlsem_vignette.Rnw:565-566
###################################################
edgeCoeff(mod0)


###################################################
### code chunk number 14: dlsem_vignette.Rnw:595-596
###################################################
causalEff(mod0,from="Job",to="Pollution",lag=seq(0,20,by=5),cumul=T)


###################################################
### code chunk number 15: dlsem_vignette.Rnw:621-622
###################################################
causalEff(mod0,from="Job",to="Pollution",cumul=T)$overall


###################################################
### code chunk number 16: dlsem_vignette.Rnw:633-635 (eval = FALSE)
###################################################
## lagPlot(mod0,path="Job*Pollution")
## lagPlot(mod0,path="Job*Consum*Pollution")


###################################################
### code chunk number 17: dlsem_vignette.Rnw:642-643 (eval = FALSE)
###################################################
## lagPlot(mod0,from="Job",to="Pollution")


