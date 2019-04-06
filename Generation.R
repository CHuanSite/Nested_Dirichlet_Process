library(plotly)
library(tidyverse)
library(MCMCpack)

####################################
#
#
#Setting for the simulation
#
#
####################################


## The weight for different distributions in the mixture
w = list()
w[[1]] = c(0.75, 0.25)
w[[2]] = c(0.55, 0.45)
w[[3]] = c(0.40, 0.30, 0.30)
w[[4]] = c(0.39, 0.29, 0.29, 0.03)

## The mean for different distirubtion in the mixture
mu = list()
mu[[1]] = c(0, 3.0)
mu[[2]] = c(0, 3.0)
mu[[3]] = c(0, -2.0, 2.0)
mu[[4]] = c(0, -2.0, 2.0, 10.0)

## The variance for different distribution in the mixture
sigma = list()
sigma[[1]] = c(1.0, 2.0)
sigma[[2]] = c(1.0, 2.0)
sigma[[3]] = c(1.0, 2.0, 2.0)
sigma[[4]] = c(1.0, 2.0, 2.0, 1.0)

## The number of samples and its nested samples
J = 20
N = 500

## The density for the four mixtures
p_density = list()
p_density[[1]] = w[[1]][1] * dnorm(seq(-10, 10, 0.1), mu[[1]][1], sqrt(sigma[[1]][1])) + 
  w[[1]][2] * dnorm(seq(-10, 10, 0.1), mu[[1]][2], sqrt(sigma[[1]][2]))
p_density[[2]] = w[[2]][1] * dnorm(seq(-10, 10, 0.1), mu[[2]][1], sqrt(sigma[[2]][1])) + 
  w[[2]][2] * dnorm(seq(-10, 10, 0.1), mu[[2]][2], sqrt(sigma[[2]][2]))
p_density[[3]] = w[[3]][1] * dnorm(seq(-10, 10, 0.1), mu[[3]][1], sqrt(sigma[[3]][1])) + 
  w[[3]][2] * dnorm(seq(-10, 10, 0.1), mu[[3]][2], sqrt(sigma[[3]][2])) + 
  w[[3]][3] * dnorm(seq(-10, 10, 0.1), mu[[3]][3], sqrt(sigma[[3]][3]))
p_density[[4]] = w[[4]][1] * dnorm(seq(-10, 10, 0.1), mu[[4]][1], sqrt(sigma[[4]][1])) + 
  w[[4]][2] * dnorm(seq(-10, 10, 0.1), mu[[4]][2], sqrt(sigma[[4]][2])) + 
  w[[4]][3] * dnorm(seq(-10, 10, 0.1), mu[[4]][3], sqrt(sigma[[4]][3])) + 
  w[[4]][4] * dnorm(seq(-10, 10, 0.1), mu[[4]][4], sqrt(sigma[[4]][4]))

plot_ly(x = seq(-10, 10, 0.1), y  = p_density[[1]], name = 'T1', type = 'scatter', mode = 'lines') %>%
  add_trace(y = p_density[[2]], name = 'T2', mode = 'lines') %>%
  add_trace(y = p_density[[3]], name = 'T3', mode = 'lines') %>%
  add_trace(y = p_density[[4]], name = 'T4', mode = 'lines') %>%
  layout(title = "True distributions used in the simulation study")

##############################################
#
#
# Generate the simulation data
#
#
###############################################


## Number of samples for each distribution
n_sample = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4)

## Do the permutation
n_sample = n_sample[sample(1 : J, replace = FALSE)]

## Store the data in dat
dat = c()

## Store the group in group
group = c()

## Store the data for each group
data_store = list()

## Generate the simulation data
for(j in 1 : J){
  dat = c()
  index_j = n_sample[j]
  print(index_j)
  ind = sample(1 : length(w[[index_j]]), N, replace = TRUE, prob = w[[index_j]])
  for(k in 1 : length(w[[index_j]])){
    dat = c(dat, rnorm(sum(ind == k), mu[[index_j]][k], sqrt(sigma[[index_j]][k])))
  }
  data_store[[j]] = dat
}

