### R code from vignette source 'dlsem_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_vignette.Rnw:375-376 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_vignette.Rnw:383-384 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_vignette.Rnw:416-417
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_vignette.Rnw:422-424
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_vignette.Rnw:465-470
###################################################
indus.code <- list(
  Job ~ 1,
  Consum~quec.lag(Job,0,15),
  Pollution~quec.lag(Job,0,15)+quec.lag(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_vignette.Rnw:526-528
###################################################
indus.global <- list(adapt=T,max.gestation=3,max.lead=15,min.width=5,sign="+")
indus.local <- list()


###################################################
### code chunk number 7: dlsem_vignette.Rnw:535-543
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
### code chunk number 8: dlsem_vignette.Rnw:548-554
###################################################
indus.global <- list(adapt=T,min.width=5)
indus.local <- list(
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 9: dlsem_vignette.Rnw:605-607
###################################################
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T)


###################################################
### code chunk number 10: dlsem_vignette.Rnw:621-622
###################################################
summary(indus.mod)


###################################################
### code chunk number 11: dlsem_vignette.Rnw:640-641 (eval = FALSE)
###################################################
## plot(indus.mod)


###################################################
### code chunk number 12: dlsem_vignette.Rnw:687-688
###################################################
causalEff(indus.mod,from="Job",to="Pollution",cumul=T)


###################################################
### code chunk number 13: dlsem_vignette.Rnw:710-712 (eval = FALSE)
###################################################
## lagPlot(indus.mod,path="Job*Pollution")
## lagPlot(indus.mod,path="Job*Consum*Pollution")


###################################################
### code chunk number 14: dlsem_vignette.Rnw:719-720 (eval = FALSE)
###################################################
## lagPlot(indus.mod,from="Job",to="Pollution")


###################################################
### code chunk number 15: dlsem_vignette.Rnw:750-768
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
### code chunk number 16: dlsem_vignette.Rnw:779-780
###################################################
lapply(list(QUEC=indus.mod,QDEC=indus.mod_2,GAMMA=indus.mod_3),BIC)


