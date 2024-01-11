library(sandwich)
library(lmtest)
library(tmvtnorm)
library(truncnorm)
library(ggplot2)
library(abind)
library(tidyverse)
rm(list = ls())

# Einstellungen
N <- 100
simulations <- 3000
B <- c(1, 1, 1, 1, 1)
n_pred <- length(B) - 1
sd = c(2,2.5,3,3.5)

B0 = 1
B1 = 1
B2 = 1
B3 = 1
B4 = 1

# Kovarianzmatrix  randomisiert erstellen
#r <- rtruncnorm(n_pred * (n_pred - 1) / 2, a = 0.2, b = 0.7, mean = 0, sd = 1)
#sigma <- diag(n_pred)
#sigma[upper.tri(sigma)] <- r
#sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
#sigma <- sigma * (outer(sd, sd))


## Matrix festhalten
sigma <- matrix(c(2,0,0,0,
                  0,2,0,0,
                  0,0,2,0,
                  0,0,0,2), nrow = 4, ncol =4 )
 
# effiziente Schätzung des Modells
function_allModels <- function(data, models = c("HC0", "HC1","HC3","HC2"), 
                               res = c("rnorm","chisq","t")) {
  
  
  predictors <- rtmvnorm(N, mean = rep(3, n_pred), sigma = 
                           sigma, lower = rep(0, n_pred), upper = rep(100, n_pred))
  
  #R^2 festhalten
  
  R2 <- 0.3
  var_yhat <- (c(B1, B2, B3, B4) %*% cov(predictors) %*% c(B1, B2, B3, B4))
  sigma_e <- as.numeric(sqrt(var_yhat*(1-R2)/R2))
  a <- sigma_e/mean(sqrt(predictors[,1]))
  
  # verschiedene Formen von Heterosked.
  #res <- rnorm(N, mean = 0, sd = a*sqrt(predictors[,1]))
  #res <- rnorm(N, mean = 0, sd = sqrt(predictors[,1]))
  #res <- rchisq(N, df = ((a*sqrt(predictors[,1]))^2)/2) - ((a*sqrt(predictors[,1]))^2/2)
  #res <- rt(N, df = 5)
  
  if ("chisq" %in% res) {
  res <- rchisq(N, df = ((a*sqrt(predictors[,1]))^2)/2) - ((a*sqrt(predictors[,1]))^2/2)
  }
  if ("rnorm" %in% res) {
    res <- rnorm(N, mean = 0, sd = a*sqrt(predictors[,1]))
  }
  
  y <- cbind(1, predictors) %*% B + res
  
  allFits <- list()
  allSummaries <- list()
  
  allFits$OLS <- lm(y ~ predictors)
  sum_fit <- summary(allFits$OLS)
  allSummaries$OLS <- sum_fit$coefficients
  
  
 
  if ("HC0" %in% models) {
    allSummaries$HC0 <- coeftest(allFits$OLS, vcov. = vcovHC(allFits$OLS, type = "HC0"))
  }

  if ("HC1" %in% models) {
    allSummaries$HC1 <- coeftest(allFits$OLS, vcov. = vcovHC(allFits$OLS, type = "HC1"))
  }
  
  if ("HC2" %in% models) {
    allSummaries$HC2 <- coeftest(allFits$OLS, vcov. = vcovHC(allFits$OLS, type = "HC2"))
  }
  
  if ("HC3" %in% models) {
    allSummaries$HC3 <- coeftest(allFits$OLS, vcov. = vcovHC(allFits$OLS, type = "HC3"))
  }
  
  
  # das wäre ein array
 summary_array <-do.call("abind", args = list(allSummaries, rev.along = 0))
  
  # Überlege, ob man die Liste oder das array ausgeben soll
  
  return(summary_array)
}

results<-sapply(1:simulations, function(i) 
  function_allModels(models=c("HC0", "HC1","HC3","HC2"),
                     res = "chisq"), simplify = "array")


coefficients_HC0 <- results[, , "HC0",]
coefficients_HC1 <- results[, , "HC1",]
coefficients_HC2 <- results[, , "HC2",]
coefficients_HC3 <- results[, , "HC3",]
coefficients_OLS <- results[, , "OLS",]
# mean von Schätzungen
results_mean_HC0  <- as.data.frame(apply(results[, , "HC0",], MARGIN = 1:2, mean))
results_mean_HC1  <- as.data.frame(apply(results[, , "HC1",], MARGIN = 1:2, mean))
results_mean_HC2  <- as.data.frame(apply(results[, , "HC2",], MARGIN = 1:2, mean))
results_mean_HC3  <- as.data.frame(apply(results[, , "HC3",], MARGIN = 1:2, mean))
results_mean_OLS  <- as.data.frame(apply(results[, , "OLS",], MARGIN = 1:2, mean))





# emp_se hinzufügen ####
results_sd_HC0  <- as.data.frame(apply(results[, , "HC0",], MARGIN = 1:2, sd))
results_sd_HC1  <- as.data.frame(apply(results[, , "HC1",], MARGIN = 1:2, sd))
results_sd_HC2  <- as.data.frame(apply(results[, , "HC2",], MARGIN = 1:2, sd))
results_sd_HC3  <- as.data.frame(apply(results[, , "HC3",], MARGIN = 1:2, sd))
results_sd_OLS  <- as.data.frame(apply(results[, , "OLS",], MARGIN = 1:2, sd))

results_mean_HC0$emp_se<-results_sd_HC0$Estimate
results_mean_HC1$emp_se<-results_sd_HC1$Estimate
results_mean_HC2$emp_se<-results_sd_HC2$Estimate
results_mean_HC3$emp_se<-results_sd_HC3$Estimate
results_mean_OLS$emp_se<-results_sd_OLS$Estimate
#### ####

#Poweranalyse

power_OLS <- coefficients_OLS[,"Pr(>|t|)",] < 0.05
results_mean_OLS$power<-rowMeans(power_OLS)
power_HC0 <- coefficients_HC0[,"Pr(>|t|)",] < 0.05
results_mean_HC0$power<-rowMeans(power_HC0)
power_HC1 <- coefficients_HC1[,"Pr(>|t|)",] < 0.05
results_mean_HC1$power<-rowMeans(power_HC1)
power_HC2 <- coefficients_HC2[,"Pr(>|t|)",] < 0.05
results_mean_HC2$power<-rowMeans(power_HC2)
power_HC3 <- coefficients_HC3[,"Pr(>|t|)",] < 0.05
results_mean_HC3$power<-rowMeans(power_HC3)





combined_data_se <- as.data.frame (rbind(
  cbind(results_mean_HC0$`Std. Error`, group = "HCO"),
  cbind(results_mean_HC1$`Std. Error`, group = "HC1"),
  cbind(results_mean_HC2$`Std. Error`, group = "HC2"),
  cbind(results_mean_HC3$`Std. Error`, group = "HC3"),
  cbind(results_mean_OLS$`Std. Error`, group = "OLS"),
  cbind(results_mean_OLS$emp_se, group = "emp_se")
))
combined_data_se$simulated_Betas<- c("B0","B1","B2","B3","B4")


combined_data_se$V1 <- as.numeric(combined_data_se$V1)


filtered_data <- combined_data_se %>%
  filter(simulated_Betas != "B0")

# Plot mit ggplot2
ggplot(filtered_data, aes(x = simulated_Betas, y = V1, color = group)) +
  geom_point() +
  geom_line(aes(group = group), linetype = "dashed") +
  labs(title = " Standardfehler in Abhängigkeit von simulierten Betas",
       x = "Simulierte Betas", y = "
       Standardfehler") +
  scale_color_manual(values = c("HCO" = "red", "HC1" = "green", "HC2" = "blue", "HC3" = "purple", "OLS" = "black",
                                "emp_se" = "orange")) +
  theme_minimal()



combined_data_power <- as.data.frame (rbind(
  cbind(results_mean_HC0$power, group = "HCO"),
  cbind(results_mean_HC1$power, group = "HC1"),
  cbind(results_mean_HC2$power, group = "HC2"),
  cbind(results_mean_HC3$power, group = "HC3"),
  cbind(results_mean_OLS$power, group = "OLS")
))
combined_data_power$simulated_Betas<- c("B0","B1","B2","B3","B4")

combined_data_power$V1 <- as.numeric(combined_data_power$V1)


filtered_data2 <- combined_data_power %>%
  filter(simulated_Betas != "B0")

# Plot mit ggplot2
ggplot(filtered_data2, aes(x = simulated_Betas, y = V1, color = group)) +
  geom_point() +
  geom_line(aes(group = group), linetype = "dashed") +
  labs(title = " Power Abhängigkeit von simulierten Betas",
       x = "Simulierte Betas", y = "
       Standardfehler") +
  scale_color_manual(values = c("HCO" = "red", "HC1" = "green", "HC2" = "blue", "HC3" = "purple", "OLS" = "black")) +
  theme_minimal()




rm(power_HC0,power_HC1,power_HC2,power_HC3,power_OLS,results_sd_HC0,results_sd_HC1,
results_sd_HC2,results_sd_HC3,results_sd_OLS,coefficients_HC0,coefficients_HC1,coefficients_HC2,
coefficients_HC3,coefficients_OLS)
























#### ####



# .3609.  47.40

((c(3,2,1,0) %*% cov(predictors) %*% c(3,2,1,0)) / 0.30) - (c(3,2,1,0) %*% cov(predictors) %*% c(3,2,1,0))
mean(30*sqrt(predictors[,1]))

R2 <- 0.3
var_yhat <- (c(3,2,1,0) %*% cov(predictors) %*% c(3,2,1,0))
sigma_e <- as.numeric(sqrt(var_yhat*(1-R2)/R2))
a <- sigma_e/mean(sqrt(predictors[,1]))


N <- 10000
predictors <- rtmvnorm(N, mean = rep(3, n_pred), sigma = sigma, lower = rep(0, n_pred), upper = rep(100, n_pred))
# BERECHNE DIE SD so, dass ein bestimmtes R^2 erreicht wird
res <- rnorm(N, mean = 0, sd = a*sqrt(predictors[,1]))
y <- cbind(1, predictors) %*% B + res

summary(lm(y ~ predictors))


# variiere beta bzw. R^2 
#   R2 = 0 und beta = 0 0 0 0  -> hier sehen wir die Alpha-Fehler-Rate
#   R2 = 0.3 0.5 0.7  beta = 1 1 1 1 -> hier sehen wir die Power

# wenn das funktioniert, dann variieren wir mal die Form der Heteroskedastizität
# also ob man sqrt oder etwas anderes nimmt oder 
# nicht-normalverteilte Residuen    


