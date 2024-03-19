

#' @title Inverse of Tukey \eqn{g}-&-\eqn{h} Transformation
#' 
#' @description
#' ..
#' 
#' @param q \link[base]{double} \link[base]{vector}, quantiles \eqn{q}
#' 
#' @param q0 \link[base]{double} \link[base]{vector}, \eqn{q_0=(q-A)/B}, for internal use to increase compute speed
#' 
#' @param A,B ..
#' 
#' @param ... ..
#' 
#' @details
#' Unfortunately, function [qGH2z], the inverse of Tukey \eqn{g}-&-\eqn{h} transformation, 
#' does not have a closed form and needs to be solved numerically.
#' 
#' For compute intensive jobs, use function [.qGH2z].
#' 
#' 
#' @returns 
#' Function [qGH2z] returns ..
#' 
#' @examples
#' z = rnorm(1e3L)
#' all.equal.numeric(.qGH2z(z2qGH(z, g = .3, h = .1), g = .3, h = .1), z)
#' all.equal.numeric(.qGH2z(z2qGH(z, g = 0, h = .1), g = 0, h = .1), z)
#' all.equal.numeric(.qGH2z(z2qGH(z, g = .2, h = 0), g = .2, h = 0), z)
#' 
#' @export
qGH2z <- function(q, q0 = (q - A)/B, A = 0, B = 1, ...) {
  # ?base::is.finite finds finite AND non-missing; as fast as `rep(TRUE, times = nq)` (where nq = length(q))
  if (!length(q0)) return(numeric()) # required by ?fitdistrplus::fitdist
  out <- q0
  qok <- is.finite(q0)
  out[qok] <- .qGH2z(q0 = q0[qok], ...)
  return(out)
}
# @note
# Inspired by `OpVaR:::gh_inv`.






