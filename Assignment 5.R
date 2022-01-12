# Assignment 5 simple hierarchical

rm(list=ls())
cat("\014")
graphics.off()

setwd("~/Desktop/bayesian_assignments/assignment_5")

# data 

#reaction times:
y = c(607, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310, 300, 382, 294, 315, 250, 320, 335, 297, 315, 316, 319, 326, 280, 330, 272, 317, 311, 336, 308, 330, 311, 304, 327, 367, 414, 411, 391, 333, 425, 416, 379, 368, 284, 260, 294, 265, 275, 270, 270, 282, 281, 283, 307, 344, 318, 326, 308, 306, 294, 342, 323, 325, 402, 219, 285, 277, 283, 271, 280, 289, 288, 199, 267, 354, 234, 270, 320, 214, 252, 234, 223, 332, 268, 286, 292, 295, 292, 275, 300, 256, 282, 319, 288, 316, 265, 306, 347, 374, 328, 305, 327, 344, 275, 218, 263, 282, 386, 307, 267, 282, 314, 328, 332, 386, 462, 368, 354, 283, 335, 264, 304, 248, 239, 288, 239, 249, 272, 289, 274, 281, 232, 279, 308, 260, 309, 312, 307, 296, 280, 331, 298, 342, 297, 297, 305, 308, 510, 490, 458, 425, 522, 927, 555, 550, 516, 548, 560, 545, 633, 496, 498, 555, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327)
J = 179 #number of total observations

#participants:
ind = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23)

data <- list(
  J=179,
  P=23,
  y=log(y), # to be able to use normal.
 # ind=ind,
  g=ind
)

library(rstan)

#model
fit <- stan(
  file="assign5_v3.stan",
  data=data,
  chains=4,
  warmup=2000,
  iter=10000,
  cores=6,
  refresh=1,
  control = list(adapt_delta = 0.8, max_treedepth = 15)
)

plot(fit)
summary("beta[,4]")

list_of_draws <- extract(fit)
print(names(list_of_draws))
hist(list_of_draws$mu[,4])

fit$"beta[4]"

# Check trace plots.
library(bayesplot)
library(stringr)
fit %>%
  mcmc_trace(
    pars = c("mu", "tau", str_c("beta[", 1:data$P, "]")),
    n_warmup = 2000,
    facet_args = list(nrow = 5, labeller = label_parsed)
  )

library(shinystan)
my_shinystan <- as.shinystan(fit)
launch_shinystan(my_shinystan)

#extract fit and make histograms
fit_ss <- extract(fit, permuted = TRUE)
fit_ss$mu
beta <- fit_ss$beta
beta_4 <- beta[, 4]
hist(exp(beta_4), breaks = 100, col = "lightblue")
abline(v = mean(exp(beta_4)), col = "blue", lwd = 2)
summary(exp(beta_4))
ci_95 <- quantile(exp(beta_4), probs = c(0.025, 0.975))
abline(v=ci_95, col ="red", lwd = 1)

# summary statistics for the dude
#mean
mean(exp(beta_4))
#mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
v <- exp(beta_4)
result <- getmode(v)
print(result)
#median
median(exp(beta_4))
#95% CI
ci_95 <- quantile(exp(beta_4), probs = c(0.025, 0.975))
print(ci_95)

# whole group mean (question 2)
mean(y)

# generated quantities for random new individual
print(fit_ss$logy_predicted)
hist(exp(fit_ss$logy_predicted), breaks = 100, col = "lightblue")
abline(v = mean(exp(fit_ss$logy_predicted)), col = "blue", lwd = 2)
ci_95 <- quantile(exp(fit_ss$logy_predicted), probs = c(0.025, 0.975))
abline(v=ci_95, col ="red", lwd = 1)

# comparing to the values on the reaction time website:
median(exp(fit_ss$logy_predicted))
mean(exp(fit_ss$logy_predicted))

# Figures comparing hierarchical model and thetas using individual sample means:

dim(fit)
dimnames(fit)
is.array(fit)
d <- as.data.frame(fit_ss)
d2<- d %>% gather("Participant", "Cases", 3:25)

# Making a figure for the last question:
# First, I want the mean log reaction time of each participant sampled.
library(dplyr)
sample <- data.frame(y, ind)
newsample <- aggregate(sample[, 1], list(sample$ind), mean)
#Now newsample$x is a list of the means for each participant.

# Next, theta values from the extracted model fit for everybody:
fit_ss$mu

smallMu <- sample(fit_ss$mu, size=1000, replace =F)

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

histSample2 <- hist(log(y), breaks=10, plot = FALSE) # Save histogram of log(y)
histThetas2 <- hist(smallMu, breaks=10, plot = FALSE) # Save histogram of a sample of mu.

plot(histSample2, col = c1, ylim=c(0,400)) # Plot 1st histogram using a transparent color
plot(histThetas2, col = c2, add = TRUE) # Add 2nd histogram using different color

