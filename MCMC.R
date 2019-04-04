##############################################################################
#
#
# Hyperparameters for the model
#
#
##############################################################################

## The number of samples and its nested samples
J = 20
N = 100

## Truncated value for the nested dirichlet process
K = 35
L = 55


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
  ind_patient_store[[1]][[j]] = sample(1 : L, N, replace = TRUE, prob = 1 : L)
}

## The list to store u,pi,v,w
## To store u
u_store = list()
u_store[[1]] = rep(0, K)
## To store pi
pi_store = list()
pi_store[[1]] = rep(0, K)
## To store v
v_store = list()
v_store[[1]] = list()
for(k in 1 : K){
  v_store[[1]][[k]] = rep(0, L)
}
## To store w
w_store = list()
w_store[[1]] = list()
for(k in 1 : K){
  w_store[[1]][[k]] = rep(0, L)
}


## Parameters for the gaussian distribution
theta_store = list()
theta_store[[1]] = list()
for(k in 1 : K){
  theta_store[[1]][k] = rnorm(L)
}

## Concentration parameters for the alpha and beta
alpha_store = 2
beta_store = 2

## Variance component for the component
psi = 1

