

#' @title Helper Functions
#' 
#' @description
#' Helper functions to be used in downstream packages.
#' 
#' @param x,q ..
#' 
#' @param q0 ..
#' 
#' @param A,B,g,h ..
#' 
#' @param log ..
#' 
#' @param interval ..
#' 
#' @param tol,maxiter ..
#' 
#' @returns
#' Returns of the helper functions are not documented, for now.
#' 
#' @name TukeyGH_helper
#' @export
.dGH <- function(x, A, B, g, h, log, interval = c(-50, 50), tol = .Machine$double.eps^.25, maxiter = 1000) {
  # not compute intensive
  # use wider `interval` since not compute intensive
  if (!(nx <- length(x))) return(double(length = 0L)) # ?fitdistrplus::fitdist will test len-0 `x`
  nA <- length(A)
  nB <- length(B)
  ng <- length(g)
  nh <- length(h)
  
  xok <- is.finite(x) # ?fitdistrplus::fitdist will test exceptions of x = c(0, 1, Inf, NaN, -1)
  
  if ((nA == 1L) && (nB == 1L) && (ng == 1L) && (nh == 1L)) {
    z <- x
    if ((h < 0) || (B < 0)) { # exception handling for ?fitdistrplus::fitdist
      z[] <- NaN
      return(z)
    }
    z[xok] <- .qGH2z(q = c(x[xok]), A = A, B = B, g = g, h = h, interval = interval, tol = tol, maxiter = maxiter)
    
  } else if ((nA == nB) && (nA == ng) && (nA == nh)) {
    #if (!all(xok)) stop('my fmx algorithm do not allow NA or Inf quantile')
    if (is.matrix(x)) {
      if (dim(x)[1L] != nA) stop('nrow of `x` do not match length of `A`')
      z <- q0 <- (x - A)/B
    } else if (is.numeric(x)) {
      z <- q0 <- tcrossprod(1/B, x) - A/B
    } else stop('illegal x: ', sQuote(class(x)[1L]))
    qok <- is.finite(q0) # not `xok` when `x` ?base::is.vector
    for (i in seq_len(nA)) {
      iok <- qok[i,]
      z[i,iok] <- .qGH2z(q0 = q0[i,iok], g = g[i], h = h[i], interval = interval, tol = tol, maxiter = maxiter)
    }
    
  } else stop('length of parameters must match')
  
  if (any(id <- is.infinite(z))) { # `z` is either vector or 'matrix'
    z[id & (z < 0)] <- interval[1L]
    z[id & (z > 0)] <- interval[2L]
  }
  
  ret_log <- -z^2/2 - log(2*pi)/2 - Deriv_z2qGH(z, B = B, g = g, h = h)
  if (log) return(ret_log)
  return(exp(ret_log))
  
}




# Derivative of [z2qGH] against `z`, on the log-scale
# inspired by ?OpVaR:::deriv_gh
# Inf in `z` *will* cause trouble
# not sure of the usage of ?base::tanh and ?base::cosh in ?gk:::Qgh_deriv
Deriv_z2qGH <- function(z, B, g, h) {
  hz2 <- h * z^2
  if (length(g) == 1L) { # length(B) == length(h) == 1L; is.vector(z, mode = 'numeric')
    if (g == 0) {
      trm2 <- 1 + hz2
    } else {
      e_gz <- exp(g*z)
      trm2 <- e_gz + h * z * (e_gz - 1)/g
    }
  } else { # length(B) == length(g) == length(h); is.matrix(z); nrow(z) = length(B)
    g1 <- (g != 0)
    z_g1 <- z[g1, , drop = FALSE]
    e_gz1 <- exp(g[g1] * z_g1)
    trm2 <- 1 + hz2 # for `g == 0`, also create 'array'
    trm2[g1,] <- e_gz1 + h[g1] * z_g1 * (e_gz1 - 1)/g[g1]
  }
  
  return(log(B) + hz2/2 + log(trm2))
}




# internal workhorse of [qGH2z]
# inverse of Tukey GH transformation; compute intensive!!!
#' @rdname TukeyGH_helper
#' @export
.qGH2z <- function(
    q, q0 = (q - A)/B, # `q` and `q0` both finite AND non-missing
    A = 0, B = 1, g = 0, h = 0, # all len-1 (`A` and `B` not needed if `q0` is provided)
    interval = c(-15, 15), # smaller `interval` for \code{QLMDe} algorithm (support of standard normal distribution)
    tol = .Machine$double.eps^.25, maxiter = 1000
) {
  
  #if (!length(q0)) return(numeric()) # required by \link[fitdistrplus]{fitdist}
  g0 <- (g == 0)
  h0 <- (h == 0)
  out <- q0
  #interval <- t.default(array(interval, dim = c(2L, length(q0))))
  
  # bound issue only in [dGH], not [qfmx]
  
  if (!g0 && !h0) { # most likely to happen in ?stats::optim; put in as the first option to save time
    out[] <- vuniroot2(y = q0 * g, f = function(z) expm1(g*z) * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    # very small `h` would cause bound-issue
    return(out)
  }
  
  if (g0 && h0) return(q0)
  
  if (!g0 && h0) { # has bound but also has explicit form!
    egz <- q0 * g + 1
    if (any(id <- (egz <= 0))) {
      out[id] <- if (g < 0) Inf else -Inf
    }
    out[!id] <- log(egz[!id]) / g
    return(out)
  }
  
  if (g0 && !h0) { # wont have the bound issue if g0
    out[] <- vuniroot2(y = q0, f = function(z) z * exp(h * z^2/2), interval = interval, tol = tol, maxiter = maxiter)
    return(out)
  }
  
}




