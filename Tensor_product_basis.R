
#' 
#' Build a tensor product basi system and plot some functions.
#' 

rm(list=ls())
cat("\014")

# Parameters -------------------------------------------------------------------

# number of basis in each domain
nbasis.1 = 5
nbasis.2 = 5

# total number of tensor product basis
nbasis.sim = nbasis.1*nbasis.2

# type of basis
basis.type = "bspline"

# 1D grid
x1.grid = seq(from=0, to=1, length=101)
x2.grid = seq(from=0, to=1, length=101)

# number of basis functions on the 2D domain
nbasis <- nbasis.1*nbasis.2

# BASIS FUNCTIONS --------------------------------------------------------------

# build 1D basis functions on each of the 1D domains
if(basis.type == "fourier"){
  basis.1 <- fda::create.fourier.basis(rangeval=range(x1.grid), nbasis=nbasis.1)
  basis.2 <- fda::create.fourier.basis(rangeval=range(x2.grid), nbasis=nbasis.2)
} else if(basis.type == "bspline"){
  basis.1 <- fda::create.bspline.basis(rangeval=range(x1.grid), nbasis=nbasis.1)
  basis.2 <- fda::create.bspline.basis(rangeval=range(x2.grid), nbasis=nbasis.2)
}

# evaluate 1D basis functions on the grid
basis.1.eval <- fda::eval.basis(evalarg=x1.grid, basisobj=basis.1)
basis.2.eval <- fda::eval.basis(evalarg=x2.grid, basisobj=basis.2)

# 2D basis functions evaluated on the grid
basis.grid.eval <- vector("list", length = nbasis)
grid.2D <- expand.grid(x1.grid, x2.grid)
counter = 1

# tensor product basis construction
for(i in 1:nbasis.1)
{
  base_i_j = matrix(nrow=length(x1.grid),ncol=length(x2.grid))
  for(j in 1:nbasis.2)
  {
    # tensor product (ma sui punti della griglia, è li stess)
    base_i_j = outer(basis.1.eval[,i],basis.2.eval[,j])
    
    basis.grid.eval[[counter]] = base_i_j
    counter = counter+1
  }
}

# Plot -------------------------------------------------------------------------

# plot some basis functions randomly choosen
indices = sample(1:nbasis, size=6)

library(rgl)
r3dDefaults$windowRect <- c(0,50, 800, 800) 
par3d(mfrow3d(2,3))

for(i in indices){
  persp3d(x1.grid, x2.grid, basis.grid.eval[[i]], col='red',
          xlab="",ylab="",zlab="")
}

