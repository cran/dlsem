### R code from vignette source 'dlsem_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_vignette.Rnw:364-365 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_vignette.Rnw:372-373 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_vignette.Rnw:405-406
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_vignette.Rnw:411-413
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_vignette.Rnw:454-459
###################################################
indus.code <- list(
  Job ~ 1,
  Consum~quec.lag(Job,0,15),
  Pollution~quec.lag(Job,0,15)+quec.lag(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_vignette.Rnw:519-521
###################################################
indus.global <- list(adapt=T,max.gestation=3,max.lead=15,min.width=5,sign="+")
indus.local <- list()


###################################################
### code chunk number 7: dlsem_vignette.Rnw:528-536
###################################################
indus.global <- list()
indus.local <- list(
  adapt=c(Consum=T,Pollution=T),
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  min.width=list(Consum=c(Job=5),Pollution=c(Job=5,Consum=5)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 8: dlsem_vignette.Rnw:541-547
###################################################
indus.global <- list(adapt=T,min.width=5)
indus.local <- list(
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 9: dlsem_vignette.Rnw:598-600
###################################################
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T)


###################################################
### code chunk number 10: dlsem_vignette.Rnw:614-615
###################################################
summary(indus.mod)


###################################################
### code chunk number 11: dlsem_vignette.Rnw:633-634 (eval = FALSE)
###################################################
## plot(indus.mod)


###################################################
### code chunk number 12: dlsem_vignette.Rnw:676-677
###################################################
causalEff(indus.mod,from="Job",to="Pollution",cumul=T)


###################################################
### code chunk number 13: dlsem_vignette.Rnw:699-701 (eval = FALSE)
###################################################
## lagPlot(indus.mod,path="Job*Pollution")
## lagPlot(indus.mod,path="Job*Consum*Pollution")


###################################################
### code chunk number 14: dlsem_vignette.Rnw:708-709 (eval = FALSE)
###################################################
## lagPlot(indus.mod,from="Job",to="Pollution")


###################################################
### code chunk number 15: dlsem_vignette.Rnw:739-757
###################################################
# model 2: quadratic decreasing lag shapes
indus.code_2 <- list(
  Job ~ 1,
  Consum~qdec.lag(Job,0,15),
  Pollution~qdec.lag(Job,0,15)+qdec.lag(Consum,0,15)
  )
indus.mod_2 <- dlsem(indus.code_2,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T)
summary(indus.mod_2)$endogenous
# model 3: gamma lag shapes
indus.code_3 <- list(
  Job ~ 1,
  Consum~gamma.lag(Job,0.5,0.5),
  Pollution~gamma.lag(Job,0.5,0.5)+gamma.lag(Consum,0.5,0.5)
  )
indus.mod_3 <- dlsem(indus.code_3,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T)
summary(indus.mod_3)$endogenous


###################################################
### code chunk number 16: dlsem_vignette.Rnw:768-769
###################################################
lapply(list(QUEC=indus.mod,QDEC=indus.mod_2,GAMMA=indus.mod_3),BIC)


