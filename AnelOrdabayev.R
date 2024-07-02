#load("Downloads/test_2406214.RData")

#PROBLEM 1
#to understand relationship between gene lenth and GC content building linear 
#regression
lreg <- lm(gc_cds ~ cds_length, data = gene_tr)
summary(lreg)$coefficients
#The coefficients are negative, showing the reverse correlation between 
#gene lenth and GC content lenth
#Thus longer genes tend to have shorter GC content
#The P-value is also small, which makes this correlation significant

#Building model for prediction GC content from gene lenth and evolutionary age
lreg_1 <- lm(gc_cds ~ cds_length + oldest_phylostratum, data = gene_tr)
summary(lreg_1)$r.squared
#About 0.01388837 (1.39%) fraction of GC variance is explained by this model
#MSE for training data
mse_tr <- mean(lreg_1$residuals^2)
mse_tr
#mse of the model is 0.007610067 (full model)
mse_length <- mean(lreg$residuals^2)
mse_length
#0.007707531 for reduced model (only gene lenth)
#Full model performs slighly better than reduced model
summary(lreg_1)
#the p-value for cds_lenght is 0.0552 which is unsignificant to prove that the 
#gene lenth is enough for good prediction of GC content

#Cross-validation on testing(new) data 
#First, predicting reduced (lreg) and full(lreg_1) models on new data
pred_gene_te <- predict(lreg, newdata = gene_te)
pred_gene_te1 <- predict(lreg_1, newdata = gene_te)

#plot testing data and predicting regressions where dark gree is indicated for
#reduced and blue for full model
plot(gene_te$gc_cds ~gene_te$cds_length)
points(gene_te$cds_length, pred_gene_te, col = "darkgreen", pch = 19)
points(gene_te$cds_length, pred_gene_te1, col = "blue", pch = 19)

#MSE on testing data
mse_te <-mean((gene_te$gc_cds - pred_gene_te)^2)
mse_te1 <-mean((gene_te$gc_cds - pred_gene_te1)^2)
mse_te
mse_te1
#mse shows the results as 0.007555024 for reduced and 0.00749084 for full model
#This results show that full model performs slightly better than reduced model
#This indicates that there is no evidence of overfitting


#PROBLEM2
#In order to predict which of the features best predict abnormal brain morphology
#the logistic regression should be used
logreg1 <- glm(gene_tr$is_ab ~ gene_tr$cds_length, family = "binomial")
logreg2 <- glm(gene_tr$is_ab ~ gene_tr$gc_cds, family = "binomial")
logreg3 <- glm(gene_tr$is_ab ~ gene_tr$oldest_phylostratum, family = "binomial")
summary(logreg1)
summary(logreg2)
summary(logreg3)

#From the summary the gene lenght and evolutional age are significant predictors for 
#abnormal brain morphology due to the low p-values, while GC content 
#cannot be a good predictor due to the higher p-value than 0.05
#as well as standard error is is much higher than the estimation of the change for an unit

#To show visually the relationship between gene age and brain morphology
boxplot(gene_tr$oldest_phylostratum ~ gene_tr$is_ab)
#From the plot the conclusion can be that the older genes don't tend to have abnormal
#brain morphology

logreg <- glm(gene_tr$is_ab ~ gene_tr$cds_length + gene_tr$gc_cds + gene_tr$oldest_phylostratum, family = "binomial")
summary(logreg)
#The gene age is the best predctor due to the low st. error and p-value, while 
#gene length also shows signficant results for prediction brain morphology
#while GC content cannot be a good independent predictor