
#'
#' Simulate a FAR(p) process and plot it.
#' 

rm(list=ls())
setwd("D:/Poli/TESI/Code/Functional-ARMA-Process")
source("Functions/simulate_FARMA.R")

# FARMA simulation -------------------------------------------------------------

## Parameters ----

# number of basis in each domain
nbasis.1.sim = 5
nbasis.2.sim = 5

# total number of tensor product basis
nbasis.sim = nbasis.1.sim*nbasis.2.sim

# type of basis
basis.type = "bspline"

# 1D grid
x1.grid = seq(from=0, to=1, length=101)
x2.grid = seq(from=0, to=1, length=101)

# number of samples
sample_size = 19
burnin = 20

# d = dimension of the underlying VARMA(1)
d = nbasis.sim

# define VAR1, VMA2, my_Sigma
VAR1=matrix(0.3,d,d); diag(VAR1)=0.8; VAR1=VAR1/norm(VAR1, type = "F")/2
VAR2=matrix(0.1,d,d); diag(VAR2)=0.5; VAR2=VAR2/norm(VAR2, type = "F")/2
VAR=array(0,c(d,d,2)); VAR[,,1]=VAR1; VAR[,,2]=VAR2

VMA1=matrix(0.3,d,d); diag(VMA1)=0.8; VMA1=VMA1/norm(VMA1, type = "F")/2
VMA2=matrix(0.1,d,d); diag(VMA2)=0.5; VMA2=VMA2/norm(VMA2, type = "F")/2
VMA=array(0,c(d,d,2)); VMA[,,1]=VMA1; VMA[,,2]=VMA2

my_Sigma=matrix(0.6,d,d)
diag(my_Sigma)=1
my_Sigma=my_Sigma/2

# set seed
simulation_seed = 1
set.seed(simulation_seed)

## simulate data ----
out <- simulate_FARMA(n = sample_size, 
                      VAR = VAR, 
                      VMA = VMA,
                      x1.grid = x1.grid, 
                      x2.grid = x2.grid, 
                      nbasis.1 = nbasis.1.sim, 
                      nbasis.2 = nbasis.2.sim, 
                      burnin = burnin,
                      sigma = my_Sigma,
                      basis.type = basis.type)

# list of matrices
Xt = out$Xt

# Plot -------------------------------------------------------------------------

## Save GIF ----

source('Functions/my_save_GIF.R')
my_save_GIF(Xt, x1.grid, x2.grid, filename="FARMA")



