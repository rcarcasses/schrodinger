#ifndef _SOLSPEC_
#define _SOLSPEC_
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <functional>
#include "common.hpp"

class SolvSpec {
  protected:
		vector<Point> potential;
		Spectrum spectrum;
		std::function<double(double)> potFunc;

	public:
		// these values can be overrided later
		double xMin  = -5;
		double xMax  = 5;
		double tol   = 1e-9;
		double h     = 0.01; //                       _
		int N        = 100;  // used only in ChebSpec  |
		double dEmin = 0.1;  // used only in Numerov    > I need to improve this
		int nPoints  = 1000; // used only in Numerov  _|

		virtual void setPotential(vector<Point>);
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
    		potential[i].x = xMin + i * (xMax - xMin) / nPoints;
    		potential[i].y = potential[i].x * potential[i].x;  // harmonic oscillator for testing
    	}
        cout << "[DEBUG] Using harmonic oscillator potential..." << endl;
    }

    nPoints = potential.size();
    xMin    = potential.front().x;
    xMax    = potential.back().x;
  	h = (xMax - xMin) / nPoints;

  	// here we use the nice feature of c++ 11, the lambda functions
  	// oh dear, this is soo cool!
  	potFunc = [this] (double x){
  		// here we have to interpolate the potential since all the testing points are not in the grid!
  		// a linear interpolation is used
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

  		double x0 = potential[index].x;
  		double x1 = potential[index + 1].x;
  		double y0 = potential[index].y;
  		double y1 = potential[index + 1].y;

  		double m = (y1 - y0) / (x1 - x0);
  		return y0 + m * (x - x0);
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
