This is an R package interface to C++ code that solves the one dimensional time independent Schrodinger equation. There are two ways to find the quantum spectrum, one using a Numerov algorithm implementation and another which uses a Chebyshev polynomial approach. Numerov routine seems to be better suited to compute the first one or two eigenvalues and eigenfunctions but for higher orders the Chebyshev routine seems to perform better (less time to get the same result). The two routines implemented provide a good way to cross check results computed with different methods.

### Dependencies
The R package is created with Rcpp. The C++ code depends on the *armadillo* awesome library, make sure you have it in your path.

### Selecting the algorithm
The default algorithm is **cheb**, which is the Chebyshev version. If you want to use the Numerov version just:

```r
setSchroMethod('numerov')
```
and so on.

The specific Chebyshev algorithm implemented is described [here](http://www.tandfonline.com/doi/full/10.1080/23311835.2015.1045223).

### Example
The workflow is like this:

```r
library(schrodinger)
setSchroMethod('cheb')  # you can ignore this since *cheb* is already the default method
setPotential(px, py)
computeSpectrum(10)
getEnergies()
getWavefunctions()
```

The arguments of the function *setPotential()* are the *x* and *y* coordinates of the potential. The numerical boundaries used by Numerov are deduced from this, so be careful, if you are for instance finding the first 100 eigenfunctions of the harmonic oscillator make sure you provide a potential range such that indeed the higher order eigenfunctions are small enough at its ends.

*computeSpectrum(N)* computes the first N eigenvalues and eigenfunctions. The routine starts by looking at the global minimum of the potential and uses it as starting point to find the zeros of an internal spectral curve definition. The scanning is controlled by the second parameter of this function, which should be always smaller than the distance between two eigenvalues.

*getEnergies()* and *getWavefunctions()* do the obvious thing.

### C++
You can use the code in *src* in your own C++ project. Make sure you are linking vs *armadillo* properly while building.
