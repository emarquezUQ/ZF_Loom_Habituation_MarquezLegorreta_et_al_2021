
#### this script is to do a binomial test on the fmr1 responses to the looms based on the wt probability of response



getwd()


##setwd("C:/Users/uqemarqu/Documents/R/behaviour analysis")
##setwd("R:/EMLPHD-Q0556/R backup 20190118/behaviour analysis")

## to import the data
Firstblock<-read.csv("loomHab_Fmr1dataset_binomialTest_forR.csv",header=T,sep=",") 


N<-colnames(Firstblock)


## looking at the structure. 
str(Firstblock)
summary(Firstblock)


### to do a binomial test
test<-binom.test(26, 42, 0.307692307692308, alternative="greater")

summary(test)

test$p.value



#### doing a loop for fmr1

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$Fmr1Escap[i], Firstblock$Fmr1N[i], Firstblock$WTp[i], alternative="greater")
  #resulttestAll[i]<-resulttest$p.value
  print(c(Firstblock$loom[i],resulttest$p.value))
}

### it seems that with this test only the 2nd (p=3.056486e-05), 3rd (p=0.03403846) 
### and 16th (p=0.03969664) loom are sig. the 11th is marginal (p=0.05599378).



#### doing a loop for hets

for (i in 1:length(Firstblock$loom)){
  
  resulttest<-binom.test(Firstblock$HetsEscap[i], Firstblock$HetsN[i], Firstblock$WTp[i], alternative="greater")
  #resulttestAll[i]<-resulttest$p.value
  print(c(Firstblock$loom[i],resulttest$p.value))
}

### it seems that with this test only the 2nd (p=0.0005916218) and 9th (p=0.0395049) loom are sig.

