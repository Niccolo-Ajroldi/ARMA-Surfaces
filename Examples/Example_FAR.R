
#'
#' Simulate a FAR(p) process and plot it.
#' 

rm(list=ls())
setwd("D:/Poli/TESI/Code/Functional-ARMA-Process")
source("Functions/simulate_FAR.R")

# FAR simulation ---------------------------------------------------------------

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

# d = dimension of the underlying VAR(1)
d=nbasis.sim 
Psi1=matrix(0.3,d,d)
diag(Psi1)=0.8
Psi1=Psi1/norm(Psi1, type = "F")*0.3
Psi=array(0,c(d,d,1))
Psi[,,1]=Psi1
my_Sigma=matrix(0.6,d,d)
diag(my_Sigma)=1
my_Sigma=my_Sigma/2

# set seed
simulation_seed = 1
set.seed(simulation_seed)

## simulate data ----
out <- simulate_FAR(n = sample_size, 
                    Psi = Psi, 
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
my_save_GIF(Xt, x1.grid, x2.grid, filename="FAR")

## Plot first 3 realizations of the process ----

color = "seagreen3"

maxx = max(sapply(Xt, max))+0.5
minn = min(sapply(Xt, min))-0.5

library(latex2exp)
x11()
par(mfrow=c(1,3))
line=-9
persp(x=x1.grid, y=x2.grid, z=Xt[[1]], col=color,
      xlab="",ylab="",zlab="",zlim=c(minn,maxx),
      ticktype='detailed')
title(TeX("$Y_{1}$"), outer = FALSE,line=line)
persp(x=x1.grid, y=x2.grid, z=Xt[[2]], col=color,
      xlab="",ylab="",zlab="",zlim=c(minn,maxx),
      ticktype='detailed')
title(TeX("$Y_{2}$"), outer = FALSE,line=line)
persp(x=x1.grid, y=x2.grid, z=Xt[[3]], col=color,
      xlab="",ylab="",zlab="",zlim=c(minn,maxx),
      ticktype='detailed')
title(TeX("$Y_{3}$"), outer = FALSE,line=line)

