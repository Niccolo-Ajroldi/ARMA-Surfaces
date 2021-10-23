
#'
#' Save a gif of a surface evolving in time. 
#' The surface should be defined on the rectangle grid.1 \times grid.2
#' 
#' Xt should be a list of matrices.
#' The index of the list refers to the time. 
#' Xt[[t]] is the surface at time t.
#' 

my_save_GIF = function(Xt, grid.1, grid.2, filename)
{
  
  library(animation)
  library(latex2exp)
  
  maxx = max(sapply(Xt, max))+0.5
  minn = min(sapply(Xt, min))-0.5
  
  n = length(Xt)
  
  # gif of the FARS
  saveGIF(
    expr = {
      pb <- progress::progress_bar$new(
        format = " Saving gif [:bar] :percent in :elapsed",
        total = n, clear = FALSE, width= 60)
      for(t in 1:n)
      {
        pb$tick()
        
        main_title = TeX(paste0("$Y_{",t,"}$"))
        persp(x=grid.1, y=grid.2, z=Xt[[t]], col="seagreen3",
              xlab="",ylab="",zlab="",zlim=c(minn,maxx),
              ticktype='detailed')
        title(main_title, outer = FALSE)
      }
    },
    movie.name = paste0(filename,".gif"),
    interval = 0.4, 
    ani.width = 1300,
    ani.height = 700
  )
}
