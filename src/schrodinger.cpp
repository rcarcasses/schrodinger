#include <Rcpp.h>
#include <string>
#include "numerov.hpp"
#include "chebspec.hpp"

using namespace Rcpp;

#define NUMEROV 1
#define CHEBYSHEV 2

SolvSpec* setSchroMethod(std::string method) {
  if(method == "numerov") {
    return new Numerov();
  }
  // use chebyshev as default method
  return new ChebSpec();
}

List getEnergiesAndIndices(SolvSpec* n) {
  vector<Mode> modes = n->getSpectrum().modes;
  NumericVector energy, index;
  for(int i = 0; i < modes.size(); i++){
    energy.push_back(modes[i].energy);
    index.push_back(modes[i].index);
  }

  return List::create(Named("energy") = energy, Named("index") = index);
}

List getPotential(SolvSpec* n) {
  vector<Point> p = n->getPotential();
  int length = p.size();
  NumericVector x;
  NumericVector y;
  for(int i = 0; i < length; i++) {
    x.push_back(p[i].x);
    y.push_back(p[i].y);
  }

  return List::create(Named("x") = x, Named("y") = y);
}

NumericVector getEnergies(SolvSpec* n) {
  vector<double> energies = n->getSpectrum().getEnergies();
  NumericVector x(energies.begin(), energies.end());

  return x;
}

List getWavefunctions(SolvSpec* n) {
  List wfs;
  vector<vector<Point>> WFs = n->getSpectrum().getWavefunctions();

  for(int i = 0; i < WFs.size(); i++) {
    int length = WFs[i].size();
    NumericVector x;
    NumericVector y;
    for(int j = 0; j < length; j++) {
      x.push_back(WFs[i][j].x);
      y.push_back(WFs[i][j].y);
    }
    List wf = List::create(Named("x") = x, Named("y") = y);
    wfs.push_back(wf);
  }

  return wfs;
}

// [[Rcpp::export]]
List computeSpectrum(NumericVector px = NumericVector(),
                     NumericVector py = NumericVector(),
                     int nEigen = 3,
                     std::string method = "cheb",
                     double dE = 0.1, double tol = 1e-9) {
  SolvSpec* n = setSchroMethod(method);
  if(px.size() != py.size()) {
    cout << "Please pass two columns with the same size for the potential" << endl;
    return List::create(Named("error") = "Please pass two columns with the same size for the potential");
  }
  // translate the potential and feed it to the specific solver
  vector<Point> potential(px.size());
  for(int i = 0; i < px.size(); i++) {
    potential[i].x = px(i);
    potential[i].y = py(i);
  }

  n->setPotential(potential);
  n->dEmin = dE;
  n->tol = tol;
  n->findSpectrum(nEigen);
  return List::create(Named("energies") = getEnergies(n),
                      Named("wfs") = getWavefunctions(n));
}
