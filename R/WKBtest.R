#' @import rootSolve

#' @export
WKB <- function(p, n) {
  v <- splinefun(p$x, p$y)
  emin <- min(p$y)
  emax <- max(p$y)
  zeroFun <- Vectorize(function(E, i) {
                Delta <- Vectorize(function(x) (E - v(x)))
                roots <- uniroot.all(Delta, c(min(p$x), max(p$x)))
                xmin <- min(roots)
                xmax <- max(roots)

                # a "wall" case
                if(length(roots) == 1)
                    xmin = min(p$x)
                if(abs(xmin - xmax) < 2 * (p$x[2] - p$x[1]))
                  return(1e3)

                # cat('x', xmin, xmax,'\n')
                integrand <- function(x) {
                  delta <- Delta(x)
                  return(sqrt(abs(delta)))
                }

                ival <- integrate(integrand, xmin, xmax)$value
                ival - pi * (i - 0.5)
  })


  unlist(lapply(1:n, function(i) max(uniroot.all(function(E) zeroFun(E, i), c(emin + 0.1, emax - 1)))))
}
