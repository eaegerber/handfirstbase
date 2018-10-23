require(MASS)

MRMSE<-matrix(0,nrow=10,ncol=5)

pd4<-read.csv("pd3modified.csv",header=T)
pd4[,5]<-as.factor(pd4[,5])
pdnew<-pd4[sample(1:458),]
i=1
data.eval<-pdnew[(1+45*(i-1)):(45*i),]
data.model<-pdnew[-((1+45*(i-1)):(45*i)),]
RMSE<-matrix(0,nrow=10,ncol=8)
for (j in 1:10){
data.test<-data.model[(41*(j-1)+1):(41*j),]
data.train<-data.model[-((1+41*(j-1)):(41*j)),]
fit1<-glm.nb(Single~offset(log(PA))+log(BMI)*Hand,data.train)
fit2<-glm.nb(Single~offset(log(PA))+log(BMI)+Hand,data.train)
fit3<-glm.nb(Single~offset(log(PA))+log(BMI)*Hand+log(RpG),data = data.train)
fit4<-glm.nb(Single~offset(log(PA))+log(BMI)+Hand+log(RpG),data.train)
fit5<-glm.nb(Single~offset(log(PA))+cBMI+Hand+log(RpG),data.train)
fit6<-glm.nb(Single~offset(log(PA))+cBMI*Hand+log(RpG),data=data.train)
fit7<-glm.nb(Single~offset(log(PA))+cBMI*Hand,data.train)
fit8<-glm.nb(Single~offset(log(PA))+cBMI+Hand,data.train)
	
rmse1<-sum((exp(predict(fit1,data.test))-data.test[,2])^2)
rmse2<-sum((exp(predict(fit2,data.test))-data.test[,2])^2)
rmse3<-sum((exp(predict(fit3,data.test))-data.test[,2])^2)
rmse4<-sum((exp(predict(fit4,data.test))-data.test[,2])^2)
rmse5<-sum((exp(predict(fit5,data.test))-data.test[,2])^2)
rmse6<-sum((exp(predict(fit6,data.test))-data.test[,2])^2)
rmse7<-sum((exp(predict(fit7,data.test))-data.test[,2])^2)
rmse8<-sum((exp(predict(fit8,data.test))-data.test[,2])^2)
	
RMSE[j,]<-c(rmse1,rmse2,rmse3,rmse4,rmse5,rmse6,rmse7,rmse8)
}
 
AIC<-c(fit1$aic,fit2$aic,fit3$aic,fit4$aic,fit5$aic,fit6$aic,fit7$aic,fit8$aic)
exp(predict(fit1,data.eval))-data.eval[,2]
    
1/2*(rank(AIC)+rank(apply(RMSE,2,mean)))

dimnames(rest)[[2]]<-c("model1","model2","model3","model4","model5","model6","model7","model8")

plot((exp(predict(fit1,data.eval))-data.eval[,2])~data.eval[,2],xlab="Number of Single in evaluation set",ylab="Residual",col="#feb24c")
plot(exp(predict(fit1,data.eval))~data.eval[,2],col="#feb24c")
plot(exp(predict(fit1,data.eval))~data.eval[,2],col="#feb24c",xlab="Number of Single in evaluation set",ylab="Predicted Single")
abline(10,0.8,col="#2ca25f")

lm(exp(predict(fit1,data.eval))~data.eval[,2])
summary(lm(exp(predict(fit1,data.eval))~data.eval[,2]))
cor(exp(predict(fit1,data.eval)),data.eval[,2])

##########################
#Deviance and Pearson TEST#
###########################
pchisq(fit1$deviance,fit1$df.residual,lower.tail=FALSE)
pchisq( sum(residuals(fit1, type="pearson")^2 ),fit1$df.residual, lower.tail=FALSE)
##########################
#Plot Deviance residual vs. link#
#Plot Pearson residual vs. link#
###########################
plot(residuals(fit2) ~predict(fit2, type="response"),xlab=expression(hat(mu)),ylab="Deviance residuals", main="(a)")
abline(h=0,col="red")
plot(residuals(fit2) ~predict(fit2, type="link"),xlab=expression(hat(eta)),ylab="Deviance residuals", main="(b)")
abline(h=0,col="red")
plot(residuals(fit2,type="response") ~predict(fit2, type="link"),xlab=expression(hat(eta)),ylab="Response residuals", main="(c)")
abline(h=0)
plot(residuals(fit2,type="response") ~predict(fit2, type="response"),xlab=expression(hat(eta)),ylab="Response residuals", main="(d)")
abline(h=0)



##########################
#Check functional form of predictors#
#########################

scatter.smooth(data$PA, data$Single,xlab="PA",ylab="Single")
scatter.smooth(data$BMI, data$Single,xlab="BMI",ylab="Single")
scatter.smooth(data$RpG,data$Single)

scatter.smooth(log(data$PA), data$Single,xlab="log(PA)",ylab="Single")
scatter.smooth(log(data$BMI), data$Single,xlab="log(BMI)",ylab="Single")
scatter.smooth(log(data$RpG), data$Single)


##########################
#Fit new model using func. predi.#
#########################
fit2<-glm.nb(Single~log(PA)+log(BMI)*Hand,data)
summary(fit2)
fit5<-glm(Single~log(PA)+log(BMI)*Hand,family=poisson,train)
summary(fit5)
fit6<-glm(Single~log(PA)*log(BMI)*Hand,family=poisson,data)
fit4<-glm.nb(Single~log(PA)+log(BMI)+Hand,data)
pchisq(fit2$deviance,fit2$df.residual,lower.tail=FALSE)
pchisq( sum(residuals(fit2, type="pearson")^2 ),fit2$df.residual, lower.tail=FALSE)
anova(fit1,fit2)

##########################
#Partial residuals for outlier#
#########################
lambda <- predict(fit2, type="response")
u <- (data$Single-lambda)/lambda +coef(fit2)[2]*log(data$PA)
scatter.smooth(log(data$PA), u,ylab="Partial residual", main="Partial residuals")

z <- predict(fit2) + (data$Single-lambda)/lambda
scatter.smooth(predict(fit2), z, ylab="Linearlized response",main="Diagnostics for link")

halfnorm( rstudent(fit2))
halfnorm( influence(fit2)$hat,col=ifelse(influence(fit2)$hat>0.06,"red","black"))
halfnorm( cooks.distance(fit2))
plot(influence(fit2)$coef[,3],col=ifelse(influence(fit2)$coef[,3]>0.15,"red","black"), xlab="Obs.No", ylab="change in MBI coef")
fit3<-glm.nb(Single~PA+log(BMI)+Hand,data)
summary(fit3)


splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/3))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

require(foreign)
require(ggplot2)
require(MASS)
#################################
# check the overdispersion#
ggplot(data, aes(Single, fill = Hand)) + geom_histogram(binwidth = 1) + 
facet_grid(Hand ~ ., margins = TRUE, scales = "free")

a<-with(data, tapply(Single, Hand, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))