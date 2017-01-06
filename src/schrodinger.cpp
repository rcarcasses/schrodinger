#include <Rcpp.h>
#include <string>
#include "numerov.hpp"
#include "chebspec.hpp"

using namespace Rcpp;

#define NUMEROV 1
#define CHEBYSHEV 2

Numerov* numerov = new Numerov();
ChebSpec* chebyshev = new ChebSpec();

// use numerov by default
SolvSpec* n = numerov;
int method = NUMEROV;

// [[Rcpp::export]]
void setSchroMethod(std::string method, int N = -1) {
  if(method == "numerov") {
    n = numerov;
    method = NUMEROV;
  } else if(method == "cheb") {
    if(N > 2)
      chebyshev->setN(N);
    n = chebyshev;
    method = CHEBYSHEV;
  }
}

// [[Rcpp::export]]
List getEnergiesAndIndices() {
  vector<Mode> modes = n->getSpectrum().modes;
  NumericVector energy, index;
  for(int i = 0; i < modes.size(); i++){
    energy.push_back(modes[i].energy);
    index.push_back(modes[i].index);
  }

  return List::create(Named("energy") = energy, Named("index") = index);
}

// [[Rcpp::export]]
void setPotential(NumericVector px = NumericVector(), NumericVector py = NumericVector()) {
  if(px.size() != py.size()) {
    cout << "Please pass two columns with the same size for the potential" << endl;
    return;
  }
  // translate the potential and feed it to numerov
  vector<Point> potential(px.size());
  for(int i = 0; i < px.size(); i++) {
    potential[i].x = px(i);
    potential[i].y = py(i);
  }

  n->setPotential(potential);
}

// [[Rcpp::export]]
List getPotential() {
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

// [[Rcpp::export]]
void computeSpectrum(int nEigen, double dE = 0.1, double tol = 1e-9) {
  n->dEmin = dE;
  n->tol = tol;
  n->findSpectrum(nEigen);
}

// [[Rcpp::export]]
NumericVector getEnergies() {
  vector<double> energies = n->getSpectrum().getEnergies();
  NumericVector x(energies.begin(), energies.end());

  return x;
}

// [[Rcpp::export]]
List getWavefunctions() {
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
