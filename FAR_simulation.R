
#'
#' Simulate a FAR(p) process and plot it.
#' 

rm(list=ls())
source("D:/Poli/TESI/Code/Time-Series-CP/FAR_2D/simulate_FAR.R")

nbasis.1.sim = 5
nbasis.2.sim = 5
nbasis.sim = nbasis.1.sim*nbasis.2.sim
basis.type = "fourier"

x1.grid = seq(from=0, to=1, length=101)
x2.grid = seq(from=0, to=1, length=101)

# Parameters
sample_size = 5

# FAR(1)  ----------------------------------------------------------------------

## Parameters ----------------

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
out <- simulate_FAR(n = 19, 
                    Psi = Psi, 
                    x1.grid = x1.grid, 
                    x2.grid = x2.grid, 
                    nbasis.1 = nbasis.1.sim, 
                    nbasis.2 = nbasis.2.sim, 
                    sigma = my_Sigma,
                    basis.type = basis.type)

# list of matrices
Xt = out$Xt

# Save GIF ----------------

source('D:/Poli/TESI/Code/Time-Series-CP/FAR_2D/my_save_GIF.R')
my_save_GIF(Xt, x1.grid, x2.grid, filename="D:/Poli/TESI/Pics/Yt")

## Plot first 4 realizations of the process ----------------

colfunc <- colorRampPalette(c("mediumseagreen", "mediumseagreen"))
colors = colfunc(3)
color = "seagreen3"

maxx = max(sapply(Xt, max))+0.5
minn = min(sapply(Xt, min))-0.5

library(latex2exp)
x11()
#name_dir = "D:/Poli/TESI/Pics/"
#png(file = paste0(name_dir,"Evolution.png"), width = 6000, height = 3000, units = "px", res = 400)
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
#dev.off()


## plot some basis functions ----------------

library(rgl)
r3dDefaults$windowRect <- c(0,50, 800, 800) 
par3d(mfrow3d(2,3))
indices = sample(1:nbasis, size=6)
for(i in indices){
   persp3d(out$grid.1, out$grid.2, out$basis.grid.eval[[i]], col='red',
           xlab="",ylab="",zlab="")
}



