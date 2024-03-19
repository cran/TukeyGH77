

#' @title Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' ..
#'  
#' @param z \link[base]{double} scalar or \link[base]{vector}, standard normal quantiles.
#' 
#' @param A,B,g,h \link[base]{double} scalar or \link[base]{vector}
#' 
#' @details
#' Function [z2qGH] transforms standard normal quantiles to Tukey's \eqn{g}-&-\eqn{h} quantiles.
#' 
#' @returns
#' Function [z2qGH] returns a \link[base]{double} scalar or \link[base]{vector}.
#' 
#' @export
z2qGH <- function(z, A = 0, B = 1, g = 0, h = 0) {
  # not compute intensive
  # z must be numeric vector (i.e. not 'matrix')
  nA <- length(A)
  nB <- length(B)
  ng <- length(g)
  nh <- length(h)
  if (nA == 1L && nB == 1L && ng == 1L && nh == 1L) return(.z2qGH(z, A = A, B = B, g = g, h = h))
  if ((nA != nB) || (nA != ng) || (nA != nh)) stop('distribution parameters must be of same length')
  g0 <- (g == 0) # vector
  h0 <- (h == 0)
  q <- exp(tcrossprod(h, z^2/2))
  q[g0,] <- tcrossprod(g0[g0], z) * q[g0, , drop = FALSE] # *not* `g[g0]`
  q[!g0,] <- expm1(tcrossprod(g[!g0], z)) / g[!g0] * q[!g0, , drop = FALSE]
  return(A + q * B)
}


# Tukey GH definition/transformation; not compute intensive
.z2qGH <- function(z, A = 0, B = 1, g = 0, h = 0) {
  # `z` must be numeric vector (i.e. not 'matrix'), all arguments `A`, `B`, `g` and `h` are len-1
  g0 <- (g == 0) # scalar
  h0 <- (h == 0)
  q <- if (g0 && h0) {
    z
  } else if (g0 && !h0) {
    z * exp(h * z^2/2) 
  } else if (!g0 && h0) {
    expm1(g*z) / g
    # if ((h == 0) && (g > 0)) then q > -1/g (when z \arrow -Inf)
    # if ((h == 0) && (g < 0)) then q < -1/g (when z \arrow Inf)
    # numerically, the threshold of epsilon (for `h`) will depend on `g`
  } else {
    expm1(g*z) / g * exp(h * z^2/2)
  }
  return(A + q * B)
}

