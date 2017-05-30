#R code to fit GLMNET model to features data and human voting, find correlations etc
# LI


# get X data (features)
featuresMatrix <-read.csv("features.csv", header=T)

# adjust dark channel so higher is better
featuresMatrix$DarkChannel = -featuresMatrix$DarkChannel
featuresMatrix$LR = -featuresMatrix$LR

# originally in the feature PercRing
# values where no edges found in blurred image were set as NaN
# will set them as 0 instead.
featuresMatrix$PercRing[is.na(featuresMatrix$PercRing)] <- 0

###########################################################################################################################
# find correlations
nImages = 100;
nMethods = 13;
avgCorrs <- matrix(0, nMethods , nImages )

for (i in 1:nImages ) {
imageIndices = (i-1)*14 + 1:14;
humanScore = featuresMatrix[imageIndices ,14]
for (j in 1:nMethods ) {
tmpScore = featuresMatrix[imageIndices ,j]
avgCorrs[j,i] = cor(tmpScore, humanScore, use="complete.obs", method="kendall");
}
}
avgCorrVec = apply(avgCorrs,1,median)
avgCorrVec[13] = median(avgCorrs[13,1:100], na.rm = TRUE)

featuresNames = names(featuresMatrix);
featuresNames = featuresNames[1:13]

par(mai=c(1.02,1.4,0.82,0.42))
barplot(avgCorrVec , col = c(rep(c(2,4), 6), 2),names.arg = featuresNames ,horiz=T, las = 1)


#spearman corr
featuresCorr<-cor(featuresMatrix, use="complete.obs", method="kendall")
humanCorrs = featuresCorr[14,1:13]
par(mai=c(1.02,1.4,0.82,0.42))
barplot(humanCorrs, col = c(rep(c(2,4), 6), 2), horiz=T, las = 1)


###########################################################################################################################
# fit model

# remove un-used columns (LR score and Human BT score)
featuresMatrix$HumanBT<- NULL
featuresMatrix$LR<- NULL


# get Y data (votes)
votes <-read.csv("votes_real_balance_all.csv", header=T)

# adjust X
nFeats <- 12;
nRows <- dim(votes)[1];
# continuous
X <- matrix( 0 ,nrow = nRows,ncol = nFeats);
# categrocial
XRanks <- matrix( 0 ,nrow = nRows,ncol = nFeats);

for (i in 1:nRows)
{
	for (j in 1:nFeats)
	{
       index_1 = (votes[i,1]-1)*14+ votes[i,2];
       index_2 = (votes[i,1]-1)*14 + votes[i,3];
       X[i,j] = featuresMatrix[index_2, j]-featuresMatrix[index_1,j];
       if (X[i,j]>0)
	      {
		XRanks[i,j] = 1;
		}
	 else {
            XRanks[i,j] = 0;
		}
 
	}
}



#library(caret)
library(glmnet)
# fit GLM model
Y <-votes[1:nRows,4];
fit1 = glmnet(X, Y, family = "binomial")
plot(fit1, xvar = "dev", label = TRUE)

# continuous model
cvfit1 = cv.glmnet(X, Y, family = "binomial", type.measure = "class", nfolds = 5, standardize = TRUE)
plot(cvfit1);
mtext("Number of Features", 3, line = 2.5)
coef(cvfit1, s = "lambda.min")
cvfit1$lambda.1se
cvfit1$lambda.min
summary(cvfit1$glmnet.fit)


# cateogrical 0 1 model - do not standardize
cvfit2= cv.glmnet(XRanks, Y, family = "binomial", type.measure = "class", nfolds = 5, standardize = FALSE)

coef(cvfit2, s = "lambda.min")
plot(cvfit2)
mtext("Number of Features", 3, line = 2.5)
cvfit2$lambda.1se
cvfit2$lambda.min


# coefficients etc
cvfit2$glmnet.fit$beta[1:12,54]

#nmin
cvfit1$glmnet.fit$beta[1:12,58]
#fairly close
cvfit1$glmnet.fit$beta[1:12,44]

#nmin
cvfit2$glmnet.fit$beta[1:12,37]
#fairly close
cvfit2$glmnet.fit$beta[1:12,25]

# cross validation error (misclassification rate)
cvfit2$cvm[37]
# standard error of the cross validation error
cvfit2$cvsd[37]
cvfit2$cvm[25]
cvfit2$cvsd[25]

# predict responses (not the 0 1 voting but the probability/response)
tmp = predict(cvfit1, newx = data.matrix(featuresMatrix), s = "lambda.min", type = "response")
featuresMatrix$NewMetric = tmp

############ find correlations etc with the predcited responses for teh new metric

nImages = 100;
nMethods = 14;
avgCorrs <- matrix(0, nMethods , nImages )

for (i in 1:nImages ) {
imageIndices = (i-1)*14 + 1:14;
humanScore = featuresMatrix[imageIndices ,14]
j=15
tmpScore = featuresMatrix[imageIndices ,j]
avgCorrs[14,i] = cor(tmpScore, humanScore, use="complete.obs", method="kendall");

}
avgCorrVec = apply(avgCorrs,1,median)
avgCorrVec[13] = median(avgCorrs[13,1:100], na.rm = TRUE)

featuresNames = names(featuresMatrix);
featuresNames = featuresNames[c(2:13,1,15)]

par(mai=c(1.02,1.4,0.82,0.42))
barplot(avgCorrVec[c(2:13,1,14)] , col = c(rep(c(2,4), 7)),names.arg = featuresNames ,horiz=T, las = 1)
cor(featuresMatrix$NewMetric.1, featuresMatrix$humanBT, use="complete.obs", method="kendall");

featuresNames = names(featuresMatrix);
featuresNames = featuresNames[c(1:13,15)]

par(mai=c(1.02,1.4,0.82,0.42))
barplot(avgCorrVec , col = c(rep(c(2,4), 7)),names.arg = featuresNames ,horiz=T, las = 1)





