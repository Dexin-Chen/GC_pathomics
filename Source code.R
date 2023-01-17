##Codes for running lasso regression
library(glmnet)
training_features <- read.table("training pathomics feature.csv", sep = ",", head = T,row.names=1)
fitx <- as.matrix(training_features[,2:189])
fity <- as.matrix(training_features[,1])
fit <- glmnet(fitx,fity,family = "binomial")
plot(fit,xvar = "lambda",label = TRUE)
abline(v=-3.2956410552,lty=3)
set.seed(1111)
cvfit <- cv.glmnet(fitx,fity,family="binomial",type.measure="deviance",nfolds=10)
plot(cvfit,xlim=c(-6.8,-2))
c(cvfit$lambda.min, cvfit$lambda.1se)
coef(cvfit,s=c(cvfit$lambda.min))

pred<-predict(cvfit, newx=fitx, type="link", s="lambda.min")
write.table(pred,file = "pred.csv",sep = ",")

#Codes for cut-point determination
library(survminer)
cutpoint <- read.table("cut-point.csv", sep = ",", head = T,row.names=1)
cut_point <- surv_cutpoint(data = cutpoint, time = "time", event = "status", variables = "Sig",
                           minprop = 0.1, progressbar = TRUE)
plot(cut_point,  legend = "")
cut_point

#Codes for Competing-risk regression
library(cmprsk)
library(aod)
training_cohort <- read.table("training cohort.csv", sep = ",", head = T,row.names=1)
fstatus <- training_cohort$Group
ftime <- training_cohort$DFS
##Competing-risk regression for PS
PS_SHR <- crr(ftime,fstatus,training_cohort$Sig)
##Competing-risk regression for Age
Age_SHR <- crr(ftime,fstatus,training_cohort$Age)
summary(Age_SHR)
##Competing-risk regression for Sex
Sex_SHR <- crr(ftime,fstatus,training_cohort$Sex)
summary(Sex_SHR)
##Competing-risk regression for CEA
CEA_SHR <- crr(ftime,fstatus,training_cohort$CEA)
summary(CEA_SHR)
##Competing-risk regression for CA 19-9
CA199_SHR<-crr(ftime,fstatus,training_cohort$CA199)
summary(CA199_SHR)
##Competing-risk regression for tumor location
Body <- (training_cohort$Location==1)+0
Antrum <- (training_cohort$Location==2)+0
Location <- cbind(Body,Antrum)
Location_SHR <- crr(ftime,fstatus,Location)
summary(Location_SHR)
wald.test(Location_SHR$var,Location_SHR$coef,Terms=1:2)
##Competing-risk regression for tumor size
Size_SHR<-crr(ftime,fstatus,training_cohort$Size)
summary(Size_SHR)
##Competing-risk regression for tumor differentiation
Differentiation_SHR <- crr(ftime,fstatus,training_cohort$Differentiation)
summary(Differentiation_SHR)
##Competing-risk regression for Lauren type
Lauren_SHR <- crr(ftime,fstatus,training_cohort$Lauren)
summary(Lauren_SHR)
##Competing-risk regression for T stage
Tstage_SHR <- crr(ftime,fstatus,training_cohort$Tstage)
summary(Tstage_SHR)
##Competing-risk regression for N stage
N1 <- (training_cohort$Nstage==1)+0
N2 <- (training_cohort$Nstage==2)+0
N3a <- (training_cohort$Nstage==3)+0
N3b <- (training_cohort$Nstage==4)+0
N <- cbind(N1,N2,N3a,N3b)
Nstage_SHR <- crr(ftime,fstatus,N)
summary(Nstage_SHR)
wald.test(Nstage_SHR$var,Nstage_SHR$coef,Terms=1:4)
##Competing-risk regression for Chemotherapy
Che_SHR <- crr(ftime,fstatus,training_cohort$Che)
summary(Che_SHR)

##Multivariate competing-risk regression
Covariate <- cbind(training_cohort$Sig,training_cohort$CA199,training_cohort$Tstage,N)
Covariate_SHR <- crr(ftime,fstatus,Covariate)
Covariate_clinic <- cbind(training_cohort$CA199,training_cohort$Tstage,N)
Covariate_clinic_SHR <- crr(ftime,fstatus,Covariate_clinic)
summary(Covariate_SHR)
summary(Covariate_clinic_SHR)
wald.test(Covariate_SHR$var,Covariate_SHR$coef,Terms=4:7) #wald test for N stage
wald.test(Covariate_clinic_SHR$var,Covariate_clinic_SHR$coef,Terms=3:6) #wald test for N stage

#Competing-risk nomogram construction
library(mstate)
training_cohort$CA199 <- factor(training_cohort$CA199)
training_cohort$Tstage <- factor(training_cohort$Tstage)
training_cohort$Nstage <- factor(training_cohort$Nstage)
Competing_nomogram <- crprep("DFS",status=training_cohort$Group,data=training_cohort,cens=0,trans=c(1,2),
                             keep=c("Sig","CA199","Tstage","Nstage"))
library(rms)
ddist<-datadist(Competing_nomogram)
options(datadist='ddist')
Model1 <- cph(Surv(Tstart,Tstop,status==1)~Sig+CA199+Tstage+Nstage,data=Competing_nomogram,
              subset=failcode==1,surv=T)
surv1 <- Survival(Model1)
nom.sur <- nomogram(Model1,fun=list(function(x)1-surv1(24,x),function(x)1-surv1(36,x),function(x)1-surv1(60,x)),
                    fun.at=list(c(0.05,0.1,seq(.2,.9,by=.1),0.95),c(0.05,0.1,seq(.2,.9,by=.1),0.95),c(0.05,0.1,seq(.2,.9,by=.1),0.95)),
                    funlabel=c("2-year Probability of peritoneal metastasis","3-year Probability of peritoneal metastasis",
                               "5-year Probability of peritoneal metastasis"),lp=F)
plot(nom.sur)

#Plotting calibration curve
library(riskRegression)
library(pec)
library(survival)
training_cohort <- read.table("training cohort.csv", sep = ",", head = T,row.names=1)
training_cohort$CA199 <- factor(training_cohort$CA199)
training_cohort$Tstage <- factor(training_cohort$Tstage)
training_cohort$Nstage <- factor(training_cohort$Nstage)
validation_cohort <- read.table("validation cohort.csv", sep = ",", head = T,row.names=1)
validation_cohort$CA199 <- factor(validation_cohort$CA199)
validation_cohort$Tstage <- factor(validation_cohort$Tstage)
validation_cohort$Nstage <- factor(validation_cohort$Nstage)
fgr1 <- FGR(Hist(DFS,Group)~Sig+CA199+Tstage+Nstage,data=training_cohort,cause=1)
fgr2 <- FGR(Hist(DFS,Group)~CA199+Tstage+Nstage,data=training_cohort,cause=1)
fgr3 <- FGR(Hist(DFS,Group)~CA199,data=training_cohort,cause=1)
fgr4 <- FGR(Hist(DFS,Group)~Tstage,data=training_cohort,cause=1)
fgr5 <- FGR(Hist(DFS,Group)~Nstage,data=training_cohort,cause=1)
fgr6 <- FGR(Hist(DFS,Group)~Sig,data=training_cohort,cause=1)
calPlot(list("Fine-Gray"=fgr1),time=,24,data=validation_cohort,legend=FALSE, col="blue",lwd=2,xlab="Predicted peritoneal metastasis probability", 
        ylab="Observed peritoneal metastasis frequencies", 
        splitMethod="BootCv",B=100,M=100,q=3,percent=F,xlim=c(0,0.8),ylim=c(0,0.8),method="quantile")
calPlot(list("Fine-Gray"=fgr1),time=,36,data=validation_cohort,legend=FALSE, col="red",lwd=2,xlab="Predicted peritoneal metastasis probability", 
        ylab="Observed peritoneal metastasis frequencies", 
        splitMethod="BootCv",B=100,M=100,q=3,percent=F,xlim=c(0,0.8),ylim=c(0,0.8),method="quantile",add = TRUE)
calPlot(list("Fine-Gray"=fgr1),time=,60,data=validation_cohort,legend=FALSE, col="green",lwd=2,xlab="Predicted peritoneal metastasis probability", 
        ylab="Observed peritoneal metastasis frequencies", 
        splitMethod="BootCv",B=100,M=100,q=3,percent=F,xlim=c(0,0.8),ylim=c(0,0.8),method="quantile",add = TRUE)
#ROC comparison
ROC_comparison_training <- Score(list("Model1"=fgr1,"Model2"=fgr2),
                        formula=Hist(DFS,Group)~1,data=training_cohort,
                        times=60,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_comparison_training)
ROC_comparison_training

ROC_comparison_validation <- Score(list("Model1"=fgr1,"Model2"=fgr2),
                                 formula=Hist(DFS,Group)~1,data=validation_cohort,
                                 times=60,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_comparison_validation)
ROC_comparison_validation

#ROC of different variables
##ROC at 2 years in training cohort
ROC_training_2 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                        formula=Hist(DFS,Group)~1,data=training_cohort,
                        times=24,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_training_2)
##ROC at 3 years in training cohort
ROC_training_3 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                        formula=Hist(DFS,Group)~1,data=training_cohort,
                        times=36,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_training_3)
##ROC at 5 years in training cohort
ROC_training_5 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                        formula=Hist(DFS,Group)~1,data=training_cohort,
                        times=60,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_training_5)
##ROC at 2 years in validation cohort
ROC_validation_2 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                        formula=Hist(DFS,Group)~1,data=validation_cohort,
                        times=24,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_validation_2)
##ROC at 3 years in validation cohort
ROC_validation_3 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                          formula=Hist(DFS,Group)~1,data=validation_cohort,
                          times=36,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_validation_3)
##ROC at 5 years in validation cohort
ROC_validation_5 <- Score(list("Model3"=fgr3,"Model4"=fgr4,"Model5"=fgr5,"Model6"=fgr6),
                          formula=Hist(DFS,Group)~1,data=validation_cohort,
                          times=60,plots="ROC",summary="risk",contrasts=TRUE,cause=1)
plotROC(ROC_validation_5)








