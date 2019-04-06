##############################################################################
#
#
# Hyperparameters for the model
#
#
##############################################################################

## The number of samples and its nested samples
J = 20
N = rep(500, J)

## Truncated value for the nested dirichlet process
K = 35
L = 55

## Number of indicators to be performed
T = 100

## Hyperparameters for the prior
## Normal-Inverse-Gamma Prior
mu_nig = 0
lambda_nig = 0.01
alpha_nig = 3
beta_nig = 1

## Gamma prior for alpha and beta
a_alpha = 1
b_alpha = 1
a_beta = 1
b_beta = 1

##############################################################################
#
#
# Initialization of the parameters
#
#
##############################################################################

## The list to store the indicator for the centers
ind_center_store = list()
ind_center_store[[1]] = sample(1 : K, J, replace = TRUE, prob = 1 : K)

## The list to store the indicator for the patients within the centers
ind_patient_store = list()
ind_patient_store[[1]] = list()
for(j in 1 : J){
  ind_patient_store[[1]][[j]] = sample(1 : L, N[j], replace = TRUE, prob = 1 : L)
}

## The list to store u,pi,v,w
## To store u
u_store = list()
u_store[[1]] = rep(0.5, K)
u_store[[1]][K] = 1

## To store pi
pi_store = list()
pi_store[[1]] = rep(0, K)
pi_store[[1]][1] = u_store[[1]][1]
for(k in 2 : K){
  pi_store[[1]][k] = u_store[[1]][k]
  for(t in 1 : (k - 1)){
    pi_store[[1]][k] = pi_store[[1]][k] * (1 - u_store[[1]][t])
  }
}

## To store v
v_store = list()
v_store[[1]] = matrix(0.5, nrow = L, ncol = K)
v_store[[1]][L, ] = 1
  
## To store w
w_store = list()
w_store[[1]] = matrix(0, nrow = L, ncol = K)
w_store[[1]][1, ] = v_store[[1]][1, ]
for(l in 2 : L){
  w_store[[1]][l, ] = v_store[[1]][l, ]
  for(t in 1 : (l - 1)){
    w_store[[1]][l, ] = w_store[[1]][l, ] * (1 - v_store[[1]][t, ])
  }
}


## Parameters for the gaussian distribution
theta_store = list()
theta_store[[1]] = matrix(0, nrow = L, ncol = K)

## Concentration parameters for the alpha and beta
alpha_store = 2
beta_store = 2

## sqrt(variance) component for the component noted as psi
psi = 1



##############################################################################
#
#
# MCMC procedure
#
#
##############################################################################

for(ite in 2 : 200){
  
  ## First compute the indicator for each medical center
  ind_center_store[[ite]] = rep(0, J)
  for(j in 1 : J){
    ## For sample j
    prod_density = rep(0, K)
    for(k in 1 : K){
      ## for possiblity k  
        for(i in 1 : length(data_store[[j]])){
          sum_density = 0
          # for item i 
          for(l in 1 : L){
            sum_density = sum_density + w_store[[ite - 1]][l, k] * 
              dnorm(data_store[[j]][i], theta_store[[ite - 1]][l,k], psi[ite - 1])
          }
          prod_density[k] = prod_density[k] + log(sum_density)
        }
        #print(prod_density[k])
        prod_density[k] = prod_density[k] + log(pi_store[[ite - 1]][k])
    }
      #print(prod_density)
      #print(j)
      prod_density = exp(prod_density - max(prod_density)) / sum(exp(prod_density - max(prod_density)))
      ind_center_store[[ite]][j] = sample(1 : K, 1, replace = FALSE, prob = prod_density)
  }
  
  ## The indicator for each group 
  ind_patient_store[[ite]] = list()
  ## For every medical center
  for(j in 1 : J){
    ## For every patient in medical center j
    index_j = ind_center_store[[ite]][j]
    ind_patient_store[[ite]][[j]] = rep(0, N[j])
    for(i in 1 : N[j]){
      nested_probability = rep(0, L)
      for(l in 1 : L){
        nested_probability[l] = log(w_store[[ite - 1]][l, index_j]) + log(dnorm(data_store[[j]][i], theta_store[[ite - 1]][l, index_j], psi[ite - 1]))
      }
      nested_probability = exp(nested_probability - max(nested_probability))
      ind_patient_store[[ite]][[j]][i] = sample(1 : L, 1, replace = FALSE, prob = nested_probability)
    }
  }
  
  ## Update the weighting u
  u_store[[ite]] = rep(0, K)
  ## Update u for the pi_k
  for(k in 1 : (K - 1)){
    m_k = sum(ind_center_store[[ite]] == k)
    m_larger_k = sum(ind_center_store[[ite]] > k)
    u_store[[ite]][k] = rbeta(1, 1 + m_k, alpha_store[ite - 1] + m_larger_k)
  }
  u_store[[ite]][K] = 1
  
  ## Update the pi based on u
  pi_store[[ite]] = rep(0, K)
  pi_store[[ite]][1] = u_store[[ite]][1]
  for(k in 2 : K){
    pi_store[[ite]][k] = u_store[[ite]][k]
    for(t in 1 : (k - 1)){
      pi_store[[ite]][k] = pi_store[[ite]][k] * (1 - u_store[[ite]][t])
    }
  }
  
  ## Update the v
  ## Initialize v_store for the index 'ite'
  v_store[[ite]] = matrix(0.5, nrow = L, ncol = K)
  v_store[[ite]][L, ] = 1
  
  ## For every Medical Center type
  for(k in 1 : K){
    ## For every Patient type
    index_k = which(ind_center_store[[ite]] == k)
    for(l in 1 : (L - 1)){
      ## For each patient in the type-k medical center
      n_lk = 0
      n_larger_lk = 0
      for(t in index_k){
        n_lk = n_lk + sum(ind_center_store[[ite]][t] == l)
        n_larger_lk = n_larger_lk + sum(ind_center_store[[ite]][t] > l)
      }
      v_store[[ite]][l, k] = rbeta(1, 1 + n_lk, beta_store[ite - 1] + n_larger_lk)
    }
  }
  
  
  ## Update w
  w_store[[ite]] = matrix(0, nrow = L, ncol = K)
  w_store[[ite]][1, ] = v_store[[ite]][1, ]
  for(l in 2 : L){
    w_store[[ite]][l, ] = v_store[[ite]][l, ]
    for(t in 1 : (l - 1)){
      w_store[[ite]][l, ] = w_store[[ite]][l, ] * (1 - v_store[[ite]][t, ])
    }
  }
  
  ## Update theta
  theta_store[[ite]] = matrix(0, nrow = L, ncol = K)
  ## For every medical center type
  for(k in 1 : K){
    ## The index of medical centers which is indexed as k
    index_k = which(ind_center_store[[ite]] == k)
    n_lk = 0
    y_lk = 0
    ## For every truncated patient type
    for(l in 1 : L){
      for(t in index_k){
        index_l = which(ind_patient_store[[ite]][[t]] == l)
        n_lk = n_lk + sum(ind_patient_store[[ite]][[t]] == l)
        y_lk = y_lk + sum(data_store[[t]][index_l])
      }
      u_temp = y_lk / (n_lk + lambda_nig)  + mu_nig * lambda_nig / (n_lk + lambda_nig)
      sigma_temp = psi[ite - 1]^2 / (n_lk + lambda_nig)
      theta_store[[ite]][l, k] = rnorm(1, u_temp, sqrt(sigma_temp))
    }
  }
  
  ## Sample the concentration parameters alpha and beta
  alpha_temp = rgamma(1, a_alpha + (K - 1), b_alpha - sum(log(1 - u_store[[ite]])[1 : (K - 1)]))
  beta_temp = rgamma(1, a_beta + K * (L - 1), b_beta - sum(log(1 - v_store[[ite]][1 : (L - 1), ])))
  
  alpha_store = c(alpha_store, alpha_temp)
  beta_store = c(beta_store, beta_temp)
  
  ## Sample the variance component for the base measure
  N_temp = 0
  for(i in 1 : length(data_store)){
      N_temp = N_temp + length(data_store[[i]])
  }
  
  y_sum_temp = 0
  for(i in 1 : length(data_store)){
    for(j in 1 : length(data_store[[i]])){
      y_sum_temp = y_sum_temp + (data_store[[i]][j] - theta_store[[ite]][ind_patient_store[[ite]][[i]][j], ind_center_store[[ite]][i]])^2
      #print(y_sum_temp)
    }
  }
  alpha_psi = N_temp / 2 + K * L / 2 + (alpha_nig + 1) * L * K - 1
  beta_psi = (y_sum_temp + 2 * K * L * beta_nig + lambda_nig * sum((theta_store[[ite]] - mu_nig) * (theta_store[[ite]] - mu_nig) )) / 2
  psi = c(psi, sqrt(rinvgamma(1, alpha_psi, beta_psi)))
  print(ite)
  print(ind_center_store[[ite]])
}



## Show the heatmap for the clustering result
cluster_matrix = matrix(0, nrow = J, ncol = J)

for(t in 1 : J){
  for(s in 1 : J){
    temp_count = 0
    for(i in 30 : 174){
      temp_count = temp_count + as.numeric(ind_center_store[[i]][t] == ind_center_store[[i]][s])
    }
    cluster_matrix[t, s] = temp_count / 14
  }
}
heatmap(1 - cluster_matrix)
