library(mice)
library(survival)
library(compareC)
library(mitools)
library(mlr)
library(rms)
library(pec)
library(dplyr)

load("data RSFs Cox models.RData")# Data are available on application but safeguarded

#Continue for Cox models, Skip to RSF for random survival forests

#Choose sample to fit models to

#Both Pre and Post menopause
df<-trainingSample[,c("bc_survtime","breastci_s","age","bmicat","oacgroup","hrtnow2","folate","wlktotim","hgtm","qcl_1","timepreg","menopinc", "ethanol" )]
dftest<-testingSample[,c("bc_survtime","breastci_s","age","bmicat","oacgroup","hrtnow2","folate","wlktotim","hgtm","qcl_1","timepreg","menopinc", "ethanol" )]

#Postmenopause
df<-trainingSample[which(trainingSample$menopinc=="postmenopause"),c("id","bc_survtime","breastci_s","age","bmicat","oacgroup","hrtnow2","folate","wlktotim","hgtm","qcl_1","timepreg", "ethanol" )]
dftest<-testingSample[which(testingSample$menopinc=="postmenopause"),c("id","bc_survtime","breastci_s","age","bmicat","oacgroup","hrtnow2","folate","wlktotim","hgtm","qcl_1","timepreg", "ethanol" )]

#set reference in weight status to be "Normal"
df$bmicat = relevel(df$bmicat, ref = "[18.5,25)")
dftest$bmicat = relevel(dftest$bmicat, ref = "[18.5,25)")

df$breastci_s[which(df$breastci_s==TRUE)]=1
df$breastci_s[which(df$breastci_s==FALSE)]=0
dftest$breastci_s[which(dftest$breastci_s==TRUE)]=1
dftest$breastci_s[which(dftest$breastci_s==FALSE)]=0


#Cox with multiple imputation

#impute data 5 times 
set.seed("123")

imputed_data<-mice(data=df,m=5,method="pmm")
test_data<-mice(data=dftest,m=5,method="pmm")

mit <- imputationList(complete(imputed_data,1),complete(imputed_data,2),complete(imputed_data,3),complete(imputed_data,4),complete(imputed_data,5))

#Choose one of the following 

#Model A: All data (menopause status included as a predictor)
modelsBmicat <- with(imputed_data,coxph(Surv(bc_survtime,breastci_s) ~ ethanol + age +bmicat+hrtnow2+folate+wlktotim+hgtm+qcl_1+timepreg+menopinc))
#Model B: Postmenopause only
modelsBmicat <- with(imputed_data,coxph(Surv(bc_survtime,breastci_s) ~ ethanol + age +bmicat+hrtnow2+folate+wlktotim+hgtm+qcl_1+timepreg))
#model C: All data with interactions between menopause status and weight status
modelsBmicat <- with(imputed_data,coxph(Surv(bc_survtime,breastci_s) ~ ethanol + age +bmicat+hrtnow2+folate+wlktotim+hgtm+qcl_1+timepreg+menopinc+bmicat*menopinc))

coef<-summary(pool(modelsBmicat),conf.int=TRUE)

#Evaluation of Cox models

##create objects to store data in
sample<-list()
testsamp<-list()
cox<-list()
pred1<-list()
Predc1<-vector("numeric")
pred2<-list()
Predc2<-vector("numeric")
IBS<-c(5)
Cind<-c(5)
cval<-matrix(nrow=100,ncol=5)
acval<-c(100)

#Cross validated C-Index and IBS

for (i in 1:5){
  #select each imputation
  sample[[i]]<-complete(imputed_data,i)
  testsamp[[i]]<-complete(test_data,i)
  #train cox model
  cox[[i]]<-cph(Surv(bc_survtime,breastci_s) ~ ethanol + age +bmicat+hrtnow2+folate+wlktotim+hgtm+qcl_1+timepreg, x=TRUE, y=TRUE, data = sample[[i]],model = TRUE,surv = TRUE)
  
  PredError <- pec(object=cox[[i]],
                   formula = Surv(bc_survtime,breastci_s) ~ ethanol + age +bmicat+hrtnow2+folate+wlktotim+hgtm+qcl_1+timepreg,
                   data=testsamp[[i]],
                   traindata = sample[[i]],
                   exact = TRUE,
                   cens.model="cox",
                   splitMethod="none",
                   B=0,
                   verbose=TRUE)
  IBS[i]=crps(PredError)[2]
  Cind[i]<-0.5*validate(cox[[i]],method="crossvalidation", B=5)[1,5]+0.5
  
  #bootstraping for C-Index boxplot
  for (j in 1:100){
    cval[j,i]<-0.5*validate(cox[[i]], B=1)[1,5]+0.5
  }
}


MeanIBS=mean(IBS)
MeanCind=mean(Cind)
acval=rowMeans(cval[,1:5])
boxplot(acval,xlab="Post- and Pre-menopause", ylab="C-Index")

#############################################################################################
################RSF below this point#########################################################

relvars<-c("bc_survtime","breastci_s", "age", "bmicat", "oacgroup", "hrtnow2","ethanol", "folate", "wlktotim", "hgtm", "qcl_1", "timepreg", "menopinc")
usedata<-UKWCS[relvars]

#split data into breast cancer incidence and no incidence
databc<-usedata[which(usedata$breastci_s == 1),]
datanbc<-usedata[which(usedata$breastci_s == 0),]

#parameters and vectors for Bootstraping and Cross-validation
B=100
CIndtr<-c(B)
CIndtest<-c(B)
CVCIndtest<-c(B)
CVCIndtr<-c(B)

#Randomly shuffle the data
usedata<-usedata[sample(nrow(usedata)),]

#Create B equally size folds
foldsbc <- cut(seq(1,nrow(databc)),breaks=B,labels=FALSE)
foldsnbc <- cut(seq(1,nrow(datanbc)),breaks=B,labels=FALSE)


#K-fold validation for C-Index
for (i in 1:B){
  testIndicesbc <- which(foldsbc==i,arr.ind=TRUE)
  testIndicesnbc <- which(foldsnbc==i,arr.ind=TRUE)
  
  testingSample<-rbind(databc[testIndicesbc,],datanbc[testIndicesnbc,])
  trainingSample<-rbind(databc[-testIndicesbc,],datanbc[-testIndicesnbc,])
  rsfTrain<-rbind(trainingSample)
  rsf<-rfsrc(Surv(bc_survtime,breastci_s)~age+bmicat+oacgroup+hrtnow2+ethanol+folate+wlktotim+hgtm+qcl_1+timepreg+menopinc, data=rsfTrain,  ntree = 50, mtry = 5, na.action = "na.impute")
  test<-predict(rsf,testingSample, na.action = "na.omit")
  CVCIndtest[i]=1-estC(test$yvar$bc_survtime,test$yvar$breastci_s,test$predicted)
  CVCIndtr[i]=1-estC(rsf$yvar$bc_survtime,rsf$yvar$breastci_s,rsf$predicted)
  
}  


#Bootstraping and box plot of C-Indices
for (i in 1:B){
  #randomly select 2/3 of each for the training sample
  indbc<-sample(floor(nrow(databc)),m*floor(nrow(databc)*2/3))#data[which(data$breastci_s == FALSE),][1:nrow(databc),]
  indnbc<-sample(floor(nrow(datanbc)),m*floor(nrow(datanbc)*2/3))
  
  #combine cases and controls in testing and training samples
  trainingSample<-rbind(databc[indbc,],datanbc[indnbc,])
  testingSample<-rbind(databc[-indbc,],datanbc[-indnbc,])
  testingSamplebc<-rbind(databc[-indbc,])
  
  
  rsfTrain<-rbind(trainingSample)
  
  rsf<-rfsrc(Surv(bc_survtime,breastci_s)~age+bmicat+oacgroup+hrtnow2+ethanol+folate+wlktotim+hgtm+qcl_1+timepreg+menopinc, data=rsfTrain,  ntree = 50, mtry = 5, na.action = "na.impute")
  
  test<-predict(rsf,testingSample, na.action = "na.omit")
  CIndtest[i]=1-estC(test$yvar$bc_survtime,test$yvar$breastci_s,test$predicted)
  CIndtr[i]=1-estC(rsf$yvar$bc_survtime,rsf$yvar$breastci_s,rsf$predicted)
}

itestingSample<-mice(data=testingSample,m=1,method="pmm")
itest<-complete(itestingSample,1)

#Brier Score
PredError <- pec(object=rsf,
                 formula = Surv(bc_survtime,breastci_s) ~ 1,
                 data=itest,
                 exact = TRUE,
                 cens.model="marginal",
                 splitMethod="cv100",
                 B=0,
                 verbose=TRUE)

crps(PredError)



##Partial dependence plots

ggvar <- gg_variable(rsf, time = 15)
pp_BMI = gg_partial_coplot(rsf, xvar = "age", groups = ggvar$bmicat, surv_type = "surv", time=15, show.plots = F)

ggplot(pp_BMI[pp_BMI$group != "[0,18.5)",], aes(x=age, y=yhat, col=group, linetype=group)) + geom_smooth(se = FALSE) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + scale_colour_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))
ggplot(pp_BMI[pp_BMI$group != "[0,18.5)",], aes(x=age, y=yhat, col=group, linetype=group)) + geom_line(se = FALSE) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + scale_colour_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))

ggsave("plot2.pdf", width = 15, height = 15, units = "cm")

pp_meno = gg_partial_coplot(rsf, xvar = "age", groups = ggvar$menopinc, surv_type = "surv", time=15, show.plots = F)

ggplot(pp_meno, aes(x=age, y=yhat, col=group, linetype=group)) + geom_smooth(se = F) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + scale_colour_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 
ggplot(pp_meno, aes(x=age, y=yhat, col=group, linetype=group)) + geom_line(se = T) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + scale_colour_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 

ggsave("plot1.pdf", width = 15, height = 15, units = "cm")

# Smoothed dependence plots

ggplot(ggvar[ggvar$bmicat != "[0,18.5)",], xvar="age", aes(x=age, y=yhat, col=bmicat)) + geom_smooth(se=F) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + scale_colour_discrete(name = "Weight status", labels = c("Normal", "Overweight", "Obese")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black"))
ggplot(ggvar, xvar="age", aes(x=age, y=yhat, col=menopinc)) + geom_smooth(se=F) + theme_bw() + labs(x = "Age (Years)", y = "Predicted survival (%) at 15 years") + scale_linetype_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + scale_colour_discrete(name = "Menopause status", labels = c("Post menopause", "Pre menopause")) + theme(legend.position = c(0.2, 0.15), legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 

ggsave("plot1_15Y_dep_plot.pdf", width = 15, height = 15, units = "cm")
ggsave("plot2_15Y_dep_plot.pdf", width = 15, height = 15, units = "cm")


#######Extracting ORs from RSF###################


###Functions returning mean odds and st. error of odds through predicted probabilities for data instances

oddsfn<-function(object,data,time){
  pred<-predict(object,data, na.action= "na.omit")
  tind<- which.min(abs(pred[['time.interest']]-time))##gives the position of the time point closest to 10
  return(mean((1-pred$survival[,tind])/pred$survival[,tind]))
  # gives predicted survival at 15 years for this partial dependence (UW) at 15 years
}

sdoddssurv<-function(object,data,time){
  pred<-predict(object,data, na.action= "na.omit")
  tind<- which.min(abs(pred[['time.interest']]-time))##gives the position of the time point closest to 10
  return(sqrt(var((1-pred$survival[,tind])/pred$survival[,tind])/nrow(pred$survival)))
}

data<-rsfTrain



###___generate all the predictions for all partial dependencies of interest__###
#these are needed to produce odds ratios 



##__PD bmi__##
##underweight
dataUW<-data
dataUW[which(dataUW$bmicat != "[0,18.5)"),]$bmicat<-"[0,18.5)"
dataUW[is.na(dataUW$bmicat),]$bmicat<-"[0,18.5)"


##normal weight
dataN<-data
dataN[which(dataN$bmicat != "[18.5,25)"),]$bmicat<-"[18.5,25)"
dataN[is.na(dataN$bmicat),]$bmicat<-"[18.5,25)"



##overweight
dataOW<-data
dataOW[which(dataOW$bmicat != "[25,30)"),]$bmicat<-"[25,30)"
dataOW[is.na(dataOW$bmicat),]$bmicat<-"[25,30)"


##obese weight
dataO<-data
dataO[which(dataO$bmicat != "[30,150)"),]$bmicat<-"[30,150)"
dataO[is.na(dataO$bmicat),]$bmicat<-"[30,150)"



##__PD bmi, menopause and age__##
##under weight pre menopause
datapreUW<-data
datapreUW[which(datapreUW$menopinc != "premenopause"),]$menopinc<-"premenopause"
datapreUW[is.na(datapreUW$menopinc),]$menopinc<-"premenopause"
datapreUW$age<-40
datapreUW[which(datapreUW$bmicat != "[0,18.5)"),]$bmicat<-"[0,18.5)"
datapreUW[is.na(datapreUW$bmicat),]$bmicat<-"[0,18.5)"



##under weight post menopause
datapostUW<-data
datapostUW[which(datapostUW$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
datapostUW[is.na(datapostUW$menopinc),]$menopinc<-"postmenopause"
datapostUW$age<-50
datapostUW[which(datapostUW$bmicat != "[0,18.5)"),]$bmicat<-"[0,18.5)"
datapostUW[is.na(datapostUW$bmicat),]$bmicat<-"[0,18.5)"


##normal pre menopause
datapreN<-data
datapreN[which(datapreN$menopinc != "premenopause"),]$menopinc<-"premenopause"
datapreN[is.na(datapreN$menopinc),]$menopinc<-"premenopause"
datapreN$age<-40
datapreN[which(datapreN$bmicat != "[18.5,25)"),]$bmicat<-"[18.5,25)"
datapreN[is.na(datapreN$bmicat),]$bmicat<-"[18.5,25)"


##normal post menopause
datapostN<-data
datapostN[which(datapostN$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
datapostN[is.na(datapostN$menopinc),]$menopinc<-"postmenopause"
datapostN$age<-50
datapostN[which(datapostN$bmicat != "[18.5,25)"),]$bmicat<-"[18.5,25)"
datapostN[is.na(datapostN$bmicat),]$bmicat<-"[18.5,25)"



##over weight pre menopause 
datapreOW<-data
datapreOW[which(datapreOW$menopinc != "premenopause"),]$menopinc<-"premenopause"
datapreOW[is.na(datapreOW$menopinc),]$menopinc<-"premenopause"
datapreOW$age <- 40
datapreOW[which(datapreOW$bmicat != "[25,30)"),]$bmicat<-"[25,30)"
datapreOW[is.na(datapreOW$bmicat),]$bmicat<-"[25,30)"

#over weight post menopause 
datapostOW<-data
datapostOW[which(datapostOW$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
datapostOW[is.na(datapostOW$menopinc),]$menopinc<-"postmenopause"
datapostOW$age<-50
datapostOW[which(datapostOW$bmicat != "[25,30)"),]$bmicat<-"[25,30)"
datapostOW[is.na(datapostOW$bmicat),]$bmicat<-"[25,30)"



##obese pre menopause
dataPreO<-data
dataPreO[which(dataPreO$menopinc != "premenopause"),]$menopinc<-"premenopause"
dataPreO[is.na(dataPreO$menopinc),]$menopinc<-"premenopause"
dataPreO$age<-40
dataPreO[which(dataPreO$bmicat != "[30,150)"),]$bmicat<-"[30,150)"
dataPreO[is.na(dataPreO$bmicat),]$bmicat<-"[30,150)"

##obese post menopause
datapostO<-data
datapostO[which(datapostO$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
datapostO[is.na(datapostO$menopinc),]$menopinc<-"postmenopause"
datapostO$age<-50
datapostO[which(datapostO$bmicat != "[30,150)"),]$bmicat<-"[30,150)"
datapostO[is.na(datapostO$bmicat),]$bmicat<-"[30,150)"



##__PD bmi and menopause__##

##under weight pre menopause
data2preUW<-data
data2preUW[which(data2preUW$menopinc != "premenopause"),]$menopinc<-"premenopause"
data2preUW[which(data2preUW$bmicat != "[0,18.5)"),]$bmicat<-"[0,18.5)"
data2preUW[is.na(data2preUW$menopinc),]$menopinc<-"premenopause"
data2preUW[is.na(data2preUW$bmicat),]$bmicat<-"[0,18.5)"




##under weight post menopause
data2postUW<-data
data2postUW[which(data2postUW$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
data2postUW[which(data2postUW$bmicat != "[0,18.5)"),]$bmicat<-"[0,18.5)"
data2postUW[is.na(data2postUW$menopinc),]$menopinc<-"postmenopause"
data2postUW[is.na(data2postUW$bmicat),]$bmicat<-"[0,18.5)"


##normal pre menopause
data2preN<-data
data2preN[which(data2preN$menopinc != "premenopause"),]$menopinc<-"premenopause"
data2preN[which(data2preN$bmicat != "[18.5,25)"),]$bmicat<-"[18.5,25)"
data2preN[is.na(data2preN$menopinc),]$menopinc<-"premenopause"
data2preN[is.na(data2preN$bmicat),]$bmicat<-"[18.5,25)"


##normal post menopause
data2postN<-data
data2postN[which(data2postN$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
data2postN[which(data2postN$bmicat != "[18.5,25)"),]$bmicat<-"[18.5,25)"
data2postN[is.na(data2postN$menopinc),]$menopinc<-"postmenopause"
data2postN[is.na(data2postN$bmicat),]$bmicat<-"[18.5,25)"

##over weight pre menopause 
data2preOW<-data
data2preOW[which(data2preOW$menopinc != "premenopause"),]$menopinc<-"premenopause"
data2preOW[which(data2preOW$bmicat != "[25,30)"),]$bmicat<-"[25,30)"
data2preOW[is.na(data2preOW$menopinc),]$menopinc<-"premenopause"
data2preOW[is.na(data2preOW$bmicat),]$bmicat<-"[25,30)"


#over weight post menopause 
data2postOW<-data
data2postOW[which(data2postOW$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
data2postOW[which(data2postOW$bmicat != "[25,30)"),]$bmicat<-"[25,30)"
data2postOW[is.na(data2postOW$menopinc),]$menopinc<-"postmenopause"
data2postOW[is.na(data2postOW$bmicat),]$bmicat<-"[25,30)"

##obese pre menopause
data2preO<-data
data2preO[which(data2preO$menopinc != "premenopause"),]$menopinc<-"premenopause"
data2preO[which(data2preO$bmicat != "[30,150)"),]$bmicat<-"[30,150)"
data2preO[is.na(data2preO$menopinc),]$menopinc<-"premenopause"
data2preO[is.na(data2preO$bmicat),]$bmicat<-"[30,150)"


##obese post menopause
data2postO<-data
data2postO[which(data2postO$menopinc != "postmenopause"),]$menopinc<-"postmenopause"
data2postO[which(data2postO$bmicat != "[30,150)"),]$bmicat<-"[30,150)"
data2postO[is.na(data2postO$menopinc),]$menopinc<-"postmenopause"
data2postO[is.na(data2postO$bmicat),]$bmicat<-"[30,150)"


##___generate odds ratios for bmi cats without age PD for table___###
odd<-matrix(nrow = 3,ncol = 4,dimnames = list(c("pre","post","total"),c("UW","N","OW","O")))
odd[1,2]<-oddsfn(rsf,data2preN,15)
odd[2,2]<-oddsfn(rsf,data2postN,15)
odd[1,3]<-oddsfn(rsf,data2preOW,15)
odd[2,3]<-oddsfn(rsf,data2postOW,15)
odd[1,1]<-oddsfn(rsf,data2preUW,15)
odd[2,1]<-oddsfn(rsf,data2postUW,15)
odd[2,4]<-oddsfn(rsf,data2postO,15)
odd[1,4]<-oddsfn(rsf,data2preO,15)
odd[3,4]<-oddsfn(rsf,dataO,15)
odd[3,1]<-oddsfn(rsf,dataUW,15)
odd[3,2]<-oddsfn(rsf,dataN,15)
odd[3,3]<-oddsfn(rsf,dataOW,15)



###___generate SDs of odds ratios for bmi cats with age PD for table___###
sodd<-matrix(nrow = 3,ncol = 4,dimnames = list(c("pre","post","total"),c("UW","N","OW","O")))
sodd[1,2]<-sdoddssurv(rsf,data2preN,15)
sodd[2,2]<-sdoddssurv(rsf,data2postN,15)
sodd[1,3]<-sdoddssurv(rsf,data2preOW,15)
sodd[2,3]<-sdoddssurv(rsf,data2postOW,15)
sodd[1,1]<-sdoddssurv(rsf,data2preUW,15)
sodd[2,1]<-sdoddssurv(rsf,data2postUW,15)
sodd[2,4]<-sdoddssurv(rsf,data2postO,15)
sodd[1,4]<-sdoddssurv(rsf,data2preO,15)
sodd[3,4]<-sdoddssurv(rsf,dataO,15)
sodd[3,1]<-sdoddssurv(rsf,dataUW,15)
sodd[3,2]<-sdoddssurv(rsf,dataN,15)
sodd[3,3]<-sdoddssurv(rsf,dataOW,15)




#generate odds ratio with normal weight status the standard
oddsratio<-odds
for (i in 1:nrow(odds)){
  ref<-odds[i,2]
  oddsratio[i,]<-odds[i,]/ref
}
list3<-list()
list3[[1]]<-oddsratio
#to get upper  confidnec intervals 
uodd<-odd-1.96*sodd

#to get lower confidnec intervals 
dodd<-odd+1.96*sodd

#generate odds ratio with normal weight status the standard
uoddr<-uodd

for (i in 1:nrow(uodd)){
  ref<-dodd[i,2]
  uoddr[i,]<-uodd[i,]/ref
}


list3[[2]]<-uoddr

#generate odds ratio with normal weight status the standard
doddr<-dodd

for (i in 1:nrow(dodd)){
  ref<-uodd[i,2]
  doddr[i,]<-dodd[i,]/ref
}


list3[[3]]<-doddr


remove (pred2preN)
remove (pred2postN)
remove (pred2preOW)
remove (pred2postOW)
remove (pred2preUW)
remove (pred2postUW)
remove (pred2preO)
remove (pred2postO)
remove (predO)
remove (predUW)
remove (predOW)
remove (predN)
remove (predpostN)
remove (predpostO)
remove (predpostOW)
remove (predpostUW)
remove (predpreN)
remove (predpreOW)
remove (predpreUW)
remove (predpreO)
remove (predPreO)
remove (data2preN)
remove (data2postN)
remove (data2preOW)
remove (data2postOW)
remove (data2preUW)
remove (data2postUW)
remove (data2preO)
remove (data2postO)
remove (dataO)
remove (dataUW)
remove (dataOW)
remove (dataN)
remove (datapostN)
remove (datapostO)
remove (datapostOW)
remove (datapostUW)
remove (datapreN)
remove (datapreOW)
remove (datapreUW)
remove (dataPreO)
