This is an R package interface to C++ code that solves the one dimensional time independent Schrodinger equation. There are two ways to find the quantum spectrum, one using a Numerov algorithm implementation and another which uses a Chebyshev polynomial approach. Numerov routine seems to be better suited to compute the first one or two eigenvalues and eigenfunctions but for higher orders the Chebyshev routine seems to perform better (less time to get the same result). The two routines implemented provide a good way to cross check results computed with different methods.

### Dependencies
The R package is created with Rcpp. The C++ code depends on the *armadillo* awesome library, make sure you have it in your path.

### Install
```r
# install devtools, ignore if already installed on your system
install.packages('devtools')
# install the library itself
devtools::install_github('rcarcasses/schrodinger')
```


The specific Chebyshev algorithm implemented is described [here](http://www.tandfonline.com/doi/full/10.1080/23311835.2015.1045223).

### Usage
An example of how to use it is the following:

```r
library(schrodinger)
# set 'cheb' method (Chebyshev) with 400 interpolation points.
chebSetN(400)
# create a potential: harmonic oscillator
x <- seq(-20, 20, len = 2000)
y <- x^2
# compute the spectrum of this potential, first 30 eigenvalues and eigenfunctions
# s$energies: eigenvalues
# s$wfs: respective eigenfunctions
s <- computeSpectrum(x, y, 30)
```

### Selecting the algorithm
The default algorithm is **cheb**, which is the Chebyshev version. If you want to use the Numerov version pass `'numerov'` as 4th parameter:

```r
s <- computeSpectrum(x, y, 10, 'numerov')
```
and so on.
### C++
You can use the code in *src* in your own C++ project. Make sure you are linking vs *armadillo* properly while building.
