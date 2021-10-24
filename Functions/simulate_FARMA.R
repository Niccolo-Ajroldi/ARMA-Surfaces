
#'
#' Functional Autoregressive Moving Average process of order p,q: FARMA(p,q)
#'
#' Function to simulate FARMA(p,q) for two-dimensional functional data
#' defined on the bidimensional grid: x1.grid \times x2.grid.
#' The stochastic process is simulated starting from a basis expansion on 
#' a tensor product basis system. It's possible to choose either a tensor 
#' product basis system of fourier basis or a tensor product system of
#' bspline basis.
#' 
#' Algorithm steps:
#' (1) build basis functions on each 1D domain
#' (2) evaluate basis functions on the 1D grids
#' (3) obtain  2D basis functions by tensor product of 1D basis functions 
#' (4) simulate coefficients of basis expansion by means of a VARMA(p) model
#' (5) build functions
#' 
#' 
#' TODO: 
#'        - add packages installation
#'        - add parameter and return description

simulate_FARMA = function(n = 100,
                          VAR = NULL,
                          VMA = NULL,
                          x1.grid=NULL,
                          x2.grid=NULL,
                          nbasis.1=10,
                          nbasis.2=10,
                          burnin = 20,
                          sigma = diag((nbasis.1*nbasis.2):1)/(nbasis.1*nbasis.2),
                          basis.type = "fourier")
{
  
  # number of basis functions on the 2D dominium (equals the dimension of the vectors)
  nbasis <- nbasis.1*nbasis.2
  
  # SIMULATE A VARMA(p,q) ------------------------------------------------------
  
  # order of the autoregressive and moving average parts (useless)
  #p=dim(VAR)[3]
  #q=dim(VMA)[3]
  
  # simulate a VARMA(p,q) for the the coefficients of basis projection
  X = multiwave::varma(N=(burnin+n), k=nbasis, VAR=VAR, VMA=VMA, cov_matrix=sigma, innov=NULL)  # N x nbasis
  X = X[(burnin+1):(burnin+n),]
  
  # BASIS FUNCTIONS ------------------------------------------------------------
  
  # if no grid is provided, I use [0,1]x[0,1]
  if(is.null(x1.grid))
    x1.grid <- seq(0, 1, by=0.05)
  if(is.null(x2.grid))
    x2.grid <- seq(0, 1, by=0.05)
  
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
  
  # considero ogni combinazione di basi
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
  
  # Function reconstruction ----------------------------------------------------
  
  # ricostruisco le funzioni a partire dai coefficienti
  
  # il processo sarà una lista con n entries
  # l'elemento t della lista sarà la funzione valutata sulla griglia al tempo t
  Xt = vector("list", length = n)
  
  pb <- progress::progress_bar$new(
    format = " Simulating FARS(p) [:bar] :percent in :elapsed",
    total = n, clear = FALSE, width= 60)
  for(t in 1:n)
  {
    pb$tick()
    Xt[[t]] = matrix(0, nrow=length(x1.grid), ncol=length(x2.grid))
    for(u in 1:length(x1.grid))
    {
      for(v in 1:length(x2.grid))
      {
        for(j in 1:nbasis)
        {
          Xt[[t]][u,v] = Xt[[t]][u,v] + basis.grid.eval[[j]][u,v] * X[t,j]
        }
      }
    }
  }
  
  
  # RETURN ---------------------------------------------------------------------
  
  my_list=list(Xt, basis.grid.eval, x1.grid, x2.grid, grid.2D)
  
  my_names=c("Xt","basis.grid.eval", "grid.1", "grid.2", "grid.2D")
  
  return(structure(.Data=my_list, names=my_names))
  
}

