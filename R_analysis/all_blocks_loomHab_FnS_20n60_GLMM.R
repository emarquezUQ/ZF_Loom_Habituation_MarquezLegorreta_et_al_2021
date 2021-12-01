


#### this is to try a GLMM again on the FnS_20n60 dataset but taking into account the fish ID

getwd()

##setwd("C:/Users/uqemarqu/Documents/E Marquez/graph_testing_post_nature_review/loom-hab_paper-progress_nov-dec_2019/behaviour analysis")
##setwd("C:/Users/Siebeck_lab.GCI-PF0EZR63/Documents/ethan lab/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("FnS_20n60_forGLMM.csv",header=T,sep=",") 


### trying the looms as ordinal. 
#Firstblock$loom<-ordered(Firstblock$loom)


##### to plot it and check the structure is correct. 
install.packages("ggplot2")
library(ggplot2)
library(lattice)

meanresp1<-aggregate(Firstblock[, 1], list(Firstblock$loom,Firstblock$cond), mean)

#qplot(meanresp1, aes(x = Group.1,y = x),xlim=c(0,10),ylim=c(0,1)) #### this one is not working yet...

##plot(meanresp1,xlim=c(0,10),ylim=c(0,1)) #### this one is not working yet...

xyplot(meanresp1$x ~ meanresp1$Group.1, groups=meanresp1$Group.2,xlim=c(0,10),ylim=c(0,1))

### OR better

ggplot(meanresp1, aes(meanresp1$Group.1, meanresp1$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp1$Group.2, group = meanresp1$Group.2))

#####

install.packages("MuMIn")  ### package to calculate the Rsquare
install.packages("lme4")  ### package to do the gl

library(lme4)
require(MuMIn)

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 


#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$speed*Firstblock$ISI
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit2)
r.squaredGLMM(fit2) 

#### comparing models. it seems that as independent variables i get a slightly higher R2m (0.1 more). but that is significant
#### I think this means that the independent (not interacting) variables model is better? 
anova(fit,fit2)




#### taking away some variables. its always better with the 3 variables. and most variance is explained by loom. as expected
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)
r.squaredGLMM(fit3) 




############### now i want to try with block2 and block3.  #########



## to import the data of the 2nd block
Firstblock<-read.csv("2nd_block_loomHab_FnS_20n60.csv",header=T,sep=",") 

meanresp1<-aggregate(Firstblock[, 1], list(Firstblock$loom,Firstblock$cond), mean)


ggplot(meanresp1, aes(meanresp1$Group.1, meanresp1$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp1$Group.2, group = meanresp1$Group.2))


### trying the looms as ordinal. 
#Firstblock$loom<-ordered(Firstblock$loom)

#install.packages("MuMIn")  ### package to calculate the Rsquare
#install.packages("lme4")  ### package to do the gl

#library(lme4)
#require(MuMIn)

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 


#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$speed*Firstblock$ISI
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit2)
r.squaredGLMM(fit2) 

#### comparing models. it seems that as independent variables i get a slightly higher R2m (0.1 more). but that is significant
#### I think this means that the independent (not interacting) variables model is better? 
anova(fit,fit2)




#### taking away some variables. its always better with the 3 variables. and most variance is explained by loom. as expected
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)
r.squaredGLMM(fit3) 



## to import the data of the 3rd block
Firstblock<-read.csv("3rd_block_loomHab_FnS_20n60.csv",header=T,sep=",") 

meanresp1<-aggregate(Firstblock[, 1], list(Firstblock$loom,Firstblock$cond), mean)


ggplot(meanresp1, aes(meanresp1$Group.1, meanresp1$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp1$Group.2, group = meanresp1$Group.2))



### trying the looms as ordinal. 
#Firstblock$loom<-ordered(Firstblock$loom)

#install.packages("MuMIn")  ### package to calculate the Rsquare
#install.packages("lme4")  ### package to do the gl

#library(lme4)
#require(MuMIn)

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)

#### calculating the Rsquare. i get marginal which epresents the variance explained by the fixed effects
#### and conditional wich is interpreted as a variance explained by the entire model,  including bothfixed and random effects
r.squaredGLMM(fit) 


#### the interaction models fails to converge... not sure why. i manage to make it converge by adding the control parameters
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$speed*Firstblock$ISI
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit2)
r.squaredGLMM(fit2) 

#### comparing models. it seems that as independent variables i get a slightly higher R2m (0.1 more). but that is significant
#### I think this means that the independent (not interacting) variables model is better? 
anova(fit,fit2)




#### taking away some variables. its always better with the 3 variables. and most variance is explained by loom. as expected
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)
r.squaredGLMM(fit3) 


