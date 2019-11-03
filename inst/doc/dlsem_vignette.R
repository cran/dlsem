### R code from vignette source 'dlsem_vignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: dlsem_vignette.Rnw:363-364 (eval = FALSE)
###################################################
## install.packages("dlsem")


###################################################
### code chunk number 2: dlsem_vignette.Rnw:371-372 (eval = FALSE)
###################################################
## update.packages("dlsem")


###################################################
### code chunk number 3: dlsem_vignette.Rnw:404-405
###################################################
require(dlsem)


###################################################
### code chunk number 4: dlsem_vignette.Rnw:410-412
###################################################
data(industry)
summary(industry)


###################################################
### code chunk number 5: dlsem_vignette.Rnw:458-463
###################################################
indus.code <- list(
  Job ~ 1,
  Consum~quec.lag(Job,0,15),
  Pollution~quec.lag(Job,0,15)+quec.lag(Consum,0,15)
  )


###################################################
### code chunk number 6: dlsem_vignette.Rnw:523-525
###################################################
indus.global <- list(adapt=T,max.gestation=3,max.lead=15,min.width=5,sign="+")
indus.local <- list()


###################################################
### code chunk number 7: dlsem_vignette.Rnw:532-540
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
### code chunk number 8: dlsem_vignette.Rnw:545-551
###################################################
indus.global <- list(adapt=T,min.width=5)
indus.local <- list(
  max.gestation=list(Consum=c(Job=3),Pollution=c(Job=3,Consum=3)),
  max.lead=list(Consum=c(Job=15),Pollution=c(Job=15,Consum=15)),
  sign=list(Consum=c(Job="+"),Pollution=c(Job="+",Consum="+"))
  )


###################################################
### code chunk number 9: dlsem_vignette.Rnw:603-605
###################################################
indus.mod <- dlsem(indus.code,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T)


###################################################
### code chunk number 10: dlsem_vignette.Rnw:620-621
###################################################
summary(indus.mod)


###################################################
### code chunk number 11: dlsem_vignette.Rnw:639-640 (eval = FALSE)
###################################################
## plot(indus.mod)


###################################################
### code chunk number 12: dlsem_vignette.Rnw:682-683
###################################################
causalEff(indus.mod,from="Job",to="Pollution",cumul=T)


###################################################
### code chunk number 13: dlsem_vignette.Rnw:705-707 (eval = FALSE)
###################################################
## lagPlot(indus.mod,path="Job*Pollution")
## lagPlot(indus.mod,path="Job*Consum*Pollution")


###################################################
### code chunk number 14: dlsem_vignette.Rnw:714-715 (eval = FALSE)
###################################################
## lagPlot(indus.mod,from="Job",to="Pollution")


###################################################
### code chunk number 15: dlsem_vignette.Rnw:745-763
###################################################
# model 2: quadratic decreasing lag shapes
indus.code_2 <- list(
  Job ~ 1,
  Consum~qdec.lag(Job,0,15),
  Pollution~qdec.lag(Job,0,15)+qdec.lag(Consum,0,15)
  )
indus.mod_2 <- dlsem(indus.code_2,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T,quiet=T)
summary(indus.mod_2)$endogenous
# model 3: gamma lag shapes
indus.code_3 <- list(
  Job ~ 1,
  Consum~gamm.lag(Job,0.5,0.5),
  Pollution~gamm.lag(Job,0.5,0.5)+gamm.lag(Consum,0.5,0.5)
  )
indus.mod_3 <- dlsem(indus.code_3,group="Region",exogenous=c("Population","GDP"),
  data=industry,global.control=indus.global,local.control=indus.local,log=T,quiet=T)
summary(indus.mod_3)$endogenous


###################################################
### code chunk number 16: dlsem_vignette.Rnw:773-774
###################################################
compareModels(list(indus.mod,indus.mod_2,indus.mod_3))


