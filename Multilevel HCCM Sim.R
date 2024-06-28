library(sandwich)
library(lmtest)
library(tmvtnorm)
library(truncnorm)
library(ggplot2)
library(abind)
library(tidyverse)
library(metRology)
library(brms)
library(posterior)
library(lme4)
library(lmerTest)
library(nlme)
library(CR2)
library(microbenchmark)
rm(list = ls())

### zweiter Versuch ####
N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80

sigma <- matrix(c(2,0,
                  0,2), nrow = 2, ncol = 2 )

sigma_u <- matrix(c(2,0,0,
                    0,2,0,
                    0,0,2), nrow = 3, ncol = 3 )



function_sim<- function() { 
  
  predictors <- rtmvnorm(N*Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  e<- rnorm(N*Clstr, mean = 0, sd = 1)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = 
                  sigma_u)
  
  x <- cbind(1, predictors)
  
  
  # Designmatrix Z als Blockmatrix erstellen
  
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(u) + e
  
  y <- as.vector(y)
  
  cluster_id <- rep(1:Clstr, each = N) #cluster ID hinzufügen
  
  # Dataframe erstellen
  sim_data <- data.frame(y = y, x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  
  # Fitte das Modell mit lmer
  fit <- lmer(y ~ x1 + x2 + (1 + x1 + x2 | cluster_id), data = sim_data)
  
  sum <- summary(fit)
  coef<- sum$coefficients
  varcorr<- VarCorr(fit)
  
  
  return(varcorr)
  
  }

results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# apply um Estimates zu mitteln
apply(results, c(1, 2), mean)

# funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)  






#### erster Versuch ####
# Settings

simulations <- 200
N <- 50
Clstr <- 80
c00 <- 1
c10 <- 1
n_pred_u <- 2

sigma <- matrix(c(2,0,
                  0,2), nrow = 2, ncol = 2 )
dat <- data.frame()
function_sim <- function() {
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_pred_u), sigma = 
                           sigma)
  
  u0 <- u[,1]
  u1 <- u[,2]
  
  x<- rnorm(N*Clstr, mean = 0, sd = 1)
  e<- rnorm(N*Clstr, mean = 0, sd = 1)
  
  B0 <- c00 + u0
  B1 <- c10 + u1
  
  B1 <- matrix(rep(B1, each = N), ncol = 1)
  B0 <- matrix(rep(B0, each = N), ncol = 1)
    
  
  # Berechnung von y 
  y <- B0 + B1 * x + e
  
  dat<- as.data.frame(cbind(y,x,e,B0,B1))
  
  cluster_id <- rep(1:Clstr, each = N)
  
  # Hinzufügen der Gruppierungsvariable zum DataFrame
  dat <- as.data.frame(cbind(y, x, e, B0, B1, cluster_id))
  
  m <- y ~  x + ( 1+ x | cluster_id)
  
  fit<-lmer(m, data=dat, REML = T)  
  summary<- summary(fit)
  varcorr<- VarCorr(fit)
 
  return(varcorr)
}

results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")


attr(results[[1]], "correlation")
attr(results[[1]], "stddev")
c(results[[1]])

allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)  

lapply(results [1:2], mean)



results$cluster_id




# nächste Aufgaben

# (2) Wir wollen flexibel eine beliebige Kombination von Level-1 und Level-2 Prädiktoren simulieren können 
# inkl. Interaktion (Tipp: model.matrix())
# (3) ist ChatGPTs Indizierungsvorschlag schneller als eine logische Abfrage nach dem Clusterindex
# (4) Wie baut man Heteroskedastizität auf Level 1 und Level 2 ein?
#

?microbenchmark


z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
z <- bdiag(z_blocks)

y <- x %*% fixed + z %*% as.vector(u) + e


microbenchmark({
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
}, times = 1)

# noch schneller
microbenchmark({z <- sparseMatrix(i = rep(1:(N * Clstr), each = ncol(x)),
                                  j = rep(1:(Clstr * ncol(x)), times = N),
                                  x = as.vector(t(x)),
                                  dims = c(N * Clstr, Clstr * ncol(x)))}, times = 1)







library(MASS)
library(lme4)
library(Matrix)
library(tmvtnorm)

N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80
N_lvl2 <- 3

sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u <- matrix(c(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() { 
  
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  lvl_2_predictors <- matrix(rnorm(Clstr * N_lvl2), ncol = N_lvl2)
  
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u*lvl_2_predictors[,1]^2)
  
  x <- cbind(1, predictors)
  
  # Beispiel für Level-2 Prädiktoren
  
  
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  sim_data <- cbind(sim_data, level_2_expanded[, -1]) # Entferne die doppelte Spalte cluster_id
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + Z1 * Z2 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(u) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(varcorr)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
apply(results, c(1, 2), mean)

# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)











#### ####

library(MASS)
library(lme4)
library(Matrix)
library(tmvtnorm)

N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80
N_lvl2 <- 3

sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u <- matrix(c(2, 0, 0,
                    0, 2, 0,
                    0, 0, 2), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() { 
  
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  u <- rtmvnorm(Clstr, mean = rep(0, n_u), sigma = sigma_u)
  
  x <- cbind(1, predictors)
  
  # Beispiel für Level-2 Prädiktoren
  level_2_predictors <- matrix(rnorm(Clstr * N_lvl2), ncol = N_lvl2)
  colnames(level_2_predictors) <- c("Z1", "Z2")
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  sim_data <- cbind(sim_data, level_2_expanded[, -1]) # Entferne die doppelte Spalte cluster_id
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + Z1 * Z2 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(u) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(varcorr)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
apply(results, c(1, 2), mean)

# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)




#### mit Heteroskedastizität ####

library(MASS)
library(lme4)
library(Matrix)
library(tmvtnorm)

N <- 50
simulations <- 200
fixed <- c(1, 1, 1)
n_pred <- length(fixed) - 1
n_u <- length(fixed)
Clstr <- 80
N_lvl2 <- 3

sigma <- matrix(c(2, 0,
                  0, 2), nrow = 2, ncol = 2)

sigma_u_base <- matrix(c(2, 0, 0,
                         0, 2, 0,
                         0, 0, 2), nrow = 3, ncol = 3)

# Funktion zur Simulation
function_sim <- function() { 
  
  predictors <- rtmvnorm(N * Clstr, mean = rep(0, n_pred), sigma = sigma)
  
  e <- rnorm(N * Clstr, mean = 0, sd = (predictors[,1])^2)
  
  # Beispiel für Level-2 Prädiktoren
  level_2_predictors <- matrix(rnorm(Clstr * N_lvl2), ncol = N_lvl2)
  colnames(level_2_predictors) <- c("Z1", "Z2", "Z3")
  
  # Skalierungsfaktor für die Varianz der zufälligen Effekte basierend auf Level-2-Prädiktoren
  scaling_factor <- exp(level_2_predictors[, 1])  # Beispiel: exponentielle Skalierung basierend auf Z1
  
  # Generiere zufällige Effekte mit skalierter Varianz
  u <- matrix(0, nrow = Clstr, ncol = n_u)
  for (i in 1:Clstr) {
    scaled_sigma_u <- sigma_u_base * scaling_factor[i]
    u[i, ] <- mvrnorm(1, mu = rep(0, n_u), Sigma = scaled_sigma_u)
  }
  
  x <- cbind(1, predictors)
  
  # Cluster-ID hinzufügen und Level-2-Prädiktoren replizieren
  cluster_id <- rep(1:Clstr, each = N)
  level_2_data <- data.frame(cluster_id = 1:Clstr, level_2_predictors)
  level_2_expanded <- level_2_data[rep(1:Clstr, each = N), ]
  
  # Kombiniere Level-1- und Level-2-Daten
  sim_data <- data.frame(y = rep(0, N * Clstr), x1 = x[, 2], x2 = x[, 3], cluster_id = cluster_id)
  sim_data <- cbind(sim_data, level_2_expanded[, -1]) # Entferne die doppelte Spalte cluster_id
  
  # Formel inklusive Interaktionen
  m <- y ~ x1 + x2 + Z1 * Z2 + (1 + x1 + x2 | cluster_id)
  
  # Designmatrix Z als Blockmatrix erstellen
  z_blocks <- lapply(1:Clstr, function(i) x[((i - 1) * N + 1):(i * N), ]) 
  z <- bdiag(z_blocks)
  
  y <- x %*% fixed + z %*% as.vector(t(u)) + e
  y <- as.vector(y)
  sim_data$y <- y
  
  # Fitte das Modell mit lmer
  fit <- lmer(m, data = sim_data)
  
  sum <- summary(fit)
  coef <- sum$coefficients
  varcorr <- VarCorr(fit)
  
  return(varcorr)
}

# Simulationen durchführen
results <- sapply(1:simulations, function(i) {
  function_sim()
}, simplify = "array")

# Funktion um Estimates zu mitteln
apply(results, c(1, 2), mean)

# Funktion um Varcorr zu mitteln
allCov <- t(sapply(results, FUN = function(iSample){
  c(iSample)
}, simplify = TRUE))

apply(allCov, MARGIN = 2, mean)

