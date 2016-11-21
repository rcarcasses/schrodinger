#ifndef _SOLSPEC_
#define _SOLSPEC_
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <functional>
#include "common.hpp"

using namespace std;
class SolvSpec {
  protected:
		vector<Point> potential;
		Spectrum spectrum;
		std::function<double(double)> potFunc;

	public:
		// these values can be overrided later
		double xMin  = -15;
		double xMax  = 15;
		double tol   = 1e-9;
		double h     = 0.01; //                       _
		int N        = 100;  // used only in ChebSpec  |
		double dEmin = 0.1;  // used only in Numerov    > I need to improve this
		int nPoints  = 1000; // used only in Numerov  _|
  	double *A = NULL, *B = NULL, *C = NULL, *D = NULL;   // spline variables

		virtual void setPotential(vector<Point>);
	  // virtual void setN(int n);// used only in ChebSpec
		vector<Point> getPotential();
		virtual void findSpectrum(int nEigen) = 0;
		Spectrum getSpectrum();
		void savePotential();
		SolvSpec() {
		  cout << "Initializing spectral solver" << endl;
		}
};

vector<Point> SolvSpec::getPotential() {
    return potential;
}

Spectrum SolvSpec::getSpectrum() {
	return spectrum;
}

void SolvSpec::setPotential(vector<Point> v = vector<Point>(1)){
	// the potential can be read at this point from a file or similar, or can be hard coded.
    if(v.size() > 1) {
        potential = v;
    } else { // show an harmonic oscillator data by default
    	potential.resize(nPoints);
    	for (int i = 0; i < nPoints; i++) {
    		potential[i].x = xMin + i * (xMax - xMin) / (nPoints - 1);
    		potential[i].y = potential[i].x * potential[i].x;  // harmonic oscillator for testing
    	}
        cout << "[DEBUG] Using harmonic oscillator potential..." << endl;
    }

    nPoints = potential.size();
    xMin    = potential.front().x;
    xMax    = potential.back().x;
  	h = (xMax - xMin) / nPoints;

	 	int N = nPoints - 1;
	  // delete previous and create new ones
	  delete[] A; A = new double[N];
	  delete[] B; B = new double[N];
	  delete[] C; C = new double[N];
	  delete[] D; D = new double[N];

	  double w[N];
	  double h[N];
	  double ftt[N+1];

	  for (int i=0; i<N; i++) {
  	    w[i] = (potential[i + 1].x-potential[i].x);
  	    h[i] = (potential[i + 1].y-potential[i].y) / w[i];
	  }

	  ftt[0] = 0;
	  for (int i=0; i<N-1; i++)
	      ftt[i+1] = 3*(h[i+1]-h[i])/(w[i+1]+w[i]);

	  ftt[N] = 0;

	  for (int i=0; i<N; i++) {
  	    A[i] = (ftt[i+1]-ftt[i])/(6*w[i]);
  	    B[i] = ftt[i]/2;
  	    C[i] = h[i]-w[i]*(ftt[i+1]+2*ftt[i])/6;
  	    D[i] = potential[i].y;
	  }

  	// here we use the nice feature of c++ 11, the lambda functions
  	// oh dear, this is soo cool!
  	potFunc = [this] (double x){
  		// here we have to interpolate the potential since all the testing points are not in the grid!
  		// first get the closest left point
  		int index = 0;
  		for (int i = 0; i < potential.size(); i++){
    			if (potential[i].x > x) {
    				index = i - 1;
    				break;
    			}
    			// if we arrived to the end then it because x is the end point.
    			index = potential.size() - 2;
  		}

  		return A[index] * pow(x - potential[index].x, 3) +
  		       B[index] * pow(x - potential[index].x, 2) +
  		       C[index] * (x - potential[index].x) +
  		       D[index];
	};
}



void SolvSpec::savePotential() {
    ofstream f("potential.dat");
	for (double x = xMin; x < xMax; x += 0.01) {
		f << x << " " << potFunc(x) << endl;
	}
	f.close();
}

#endif
