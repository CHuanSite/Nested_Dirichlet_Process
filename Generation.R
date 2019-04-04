## The weight for different distributions in the mixture
w = list()
w[[1]] = c(0.75, 0.25)
w[[2]] = c(0.55, 0.45)
w[[3]] = c(0.40, 0.30, 0.30)
w[[4]] = c(0.39, 0.39, 0.29, 0.03)

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
N = 100

## The density for the four mixtures
p_density = list()
p_density[[1]] = w[[1]][1] * dnorm(seq(-10, 10, 0.1), mu[[1]][1], sqrt(sigma[[1]][1])) + 
  w[[1]][2] * dnorm(seq(-10, 10, 0.1), mu[[1]][2], sqrt(sigma[[1]][2]))

