###################################################
## perlin_noise function 
##
## Source: http://stackoverflow.com/questions/15387328/realistic-simulated-elevation-data-in-r-perlin-noise
## Tutorial: http://www.redblobgames.com/articles/noise/2d/
###################################################
perlin_noise <- function( 
  n = 5,   m = 7,    # Size of the grid for the vector field
  N = 100, M = 100   # Dimension of the image
) {
  # For each point on this n*m grid, choose a unit 1 vector
  vector_field <- apply(
    array( rnorm( 2 * n * m ), dim = c(2,n,m) ),
    # array( rexp( 2 * n * m ), dim = c(2,n,m) ),
    2:3,
    function(u) u / sqrt(sum(u^2))
  )
  f <- function(x,y) {
    # Find the grid cell in which the point (x,y) is
    i <- floor(x)
    j <- floor(y)
    stopifnot( i >= 1 || j >= 1 || i < n || j < m )
    # The 4 vectors, from the vector field, at the vertices of the square
    v1 <- vector_field[,i,j]
    v2 <- vector_field[,i+1,j]
    v3 <- vector_field[,i,j+1]
    v4 <- vector_field[,i+1,j+1]
    # Vectors from the point to the vertices
    u1 <- c(x,y) - c(i,j)
    u2 <- c(x,y) - c(i+1,j)
    u3 <- c(x,y) - c(i,j+1)
    u4 <- c(x,y) - c(i+1,j+1)
    # Scalar products
    a1 <- sum( v1 * u1 )
    a2 <- sum( v2 * u2 )
    a3 <- sum( v3 * u3 )
    a4 <- sum( v4 * u4 )
    # Weighted average of the scalar products
    s <- function(p) 3 * p^2 - 2 * p^3
    p <- s( x - i )
    q <- s( y - j )
    b1 <- (1-p)*a1 + p*a2
    b2 <- (1-p)*a3 + p*a4
    (1-q) * b1 + q * b2
  }
  xs <- seq(from = 1, to = n, length = N+1)[-(N+1)]
  ys <- seq(from = 1, to = m, length = M+1)[-(M+1)]
  outer( xs, ys, Vectorize(f) )
}


######################################
## Function to generate and normalize a perlin noise map
##
## @param n Perlin noise parameter n: N and M together define the smoothness and number of local maxes on the landscape
## @param m Perlin noise parameter m: N and M together define the smoothness and number of local maxes on the landscape
## @param dimx: x-dimension of map
## @param dimy: y-dimension of map
## @param exponent factor: how extrem do you want the global max to be?
## @param seed : Random seed to use
## @param plot: create 3D plot of map as well?
## @return matrix
######################################
getMap <- function(n=2, m=2, dimx=511, dimy=511, exponent=1, seed=1, plot=FALSE) {
  set.seed(seed)
  map <- perlin_noise(n, m, dimx, dimy)
  map <- map + abs(min(map))		# Make sure there are no negative values in the map
  map <- map^exponent				# Amplify peaks and valeys through exponential scaling
  map <- ((map - min(map)) / diff(range(map)) ) 		# Re-scale to unit interval
  
  if(plot) {
    nbcol <- 100
    color <- rev(rainbow(nbcol, start=0/6, end=4/6))
    zcol  <- cut(map, nbcol)
    persp3d(1:dim(map)[1], 1:dim(map)[2], map, col=color[zcol], xlab="", ylab="", zlab="")
  }
  return(map)
}