#### Emmanuel Marquez-Legorreta et al, 2021 
#### Brain-wide visual habituation networks in wild type and fmr1 zebrafish 
#### Nature Communications. 

#### This script is to do a GLMM on the free-swimming behavioural dataset of loom habituation changing speed and ISI. 

getwd()

setwd("R:/LOOMHAB1-Q0555/MarquezLegorreta2021Datasets/Speed_n_ISI_dataset/free-swimming_behavior/data_and_scripts")

## to import the data
Firstblock<-read.csv("FnS_20n60_forGLMM.csv",header=T,sep=",") 


### trying the looms as ordinal. 
#Firstblock$loom<-ordered(Firstblock$loom)


##### to plot it and check the structure is correct. 
install.packages("ggplot2")
library(ggplot2)
library(lattice)

meanresp1<-aggregate(Firstblock[, 1], list(Firstblock$loom,Firstblock$cond), mean)

xyplot(meanresp1$x ~ meanresp1$Group.1, groups=meanresp1$Group.2,xlim=c(0,10),ylim=c(0,1))

### OR better

ggplot(meanresp1, aes(meanresp1$Group.1, meanresp1$x),xlim=c(0,10),ylim=c(0,1)) +
  geom_line(aes(linetype = meanresp1$Group.2, group = meanresp1$Group.2))

#####

install.packages("MuMIn")  ### package to calculate the Rsquare
install.packages("lme4")  ### package to do the glmm

library(lme4)
require(MuMIn)

fit <- glmer(Firstblock$resp~Firstblock$loom+Firstblock$speed+Firstblock$ISI
             +  (1|fishID),
             family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
             data=Firstblock)


summary(fit)
# 
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: Firstblock$resp ~ Firstblock$loom + Firstblock$speed + Firstblock$ISI +  
#   (1 | fishID)
# Data: Firstblock
# Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))
# 
# AIC      BIC   logLik deviance df.resid 
# 1419.4   1445.8   -704.7   1409.4     1435 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.4399 -0.4677  0.3268  0.5160  2.1694 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# fishID (Intercept) 0.9232   0.9608  
# Number of obs: 1440, groups:  fishID, 144
# 
# Fixed effects:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           3.24765    0.26708  12.160  < 2e-16 ***
#   Firstblock$loom      -0.25365    0.02616  -9.696  < 2e-16 ***
#   Firstblock$speedslow -1.23839    0.22139  -5.594 2.22e-08 ***
#   Firstblock$ISI60ISI   0.45089    0.21737   2.074    0.038 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#   (Intr) Frstblck$l Frstblck$s
# Frstblck$lm -0.676                      
# Frstblck$sp -0.536  0.105               
# Fr$ISI60ISI -0.343 -0.041     -0.027 
# 





#### calculating the Rsquare. Marginal  represents the variance explained by the fixed effects
#### and conditional is interpreted as a variance explained by the entire model,  including both fixed and random effects
r.squaredGLMM(fit) 

#                 R2m       R2c
# theoretical 0.1864698 0.3647348
# delta       0.1347685 0.2636072


##### testing and comparing other possible models ###


#### Testing interactive model
fit2 <- glmer(Firstblock$resp~Firstblock$loom*Firstblock$speed*Firstblock$ISI
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit2)
r.squaredGLMM(fit2) 

#### comparing models. it seems that as independent variables gets a slightly higher R2m (0.1 more). but that is significantly different from the interactive model
anova(fit,fit2)



#### taking away some variables. its always better with the 3 variables. and most variance is explained by loom.
fit3 <- glmer(Firstblock$resp~Firstblock$ISI+Firstblock$speed
              +  (1|fishID),
              family="binomial", 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
              data=Firstblock)


summary(fit3)
r.squaredGLMM(fit3) 


