#ifndef _NUMEROV_
#define _NUMEROV_
#include "solvspec.hpp"

class Numerov : public SolvSpec {
	private:
		vector<Range> zeros;
		vector<Point> solLR;
		vector<Point> solRL;
		vector<Point> sol;
		function<double(double)> diffFunc = [this](double E){
			return diff(E);
		};
		Point minPot;
		void findMinPot();
		void buildSol() ;
		int getSolIndex();
		double diff(double E);
		//void plotXY(std::vector<double> x, std::vector<double> y);
		void scanForZeroRanges(int nZeros) ;

	public:
		virtual void findSpectrum(int nEigen);
};

// just returns the overall minima of the potential
void Numerov::findMinPot(){
	minPot = potential[0];

	for (int i = 0; i < potential.size(); i++)
		if (minPot.y > potential[i].y)
			minPot = potential[i];

    // cout << "minPlot " << minPot.y << " at " << minPot.x <<endl;
}

double Numerov::diff(double E){
	  // first we find the matching point
    function<double(double)> shiftedPot = [this, &E](double x) {
  		//cout<<"E=" << E << " shiftedPot(" << x <<")=" << E - potFunc(x) << endl;
  		return E - SolvSpec::potFunc(x);
  	};

  	// since we want to find preferently the right turning point we
  	// look from the minimum position to the right
  	double matchPoint = zbrent(shiftedPot, minPot.x, xMax, tol, true);
  	//cout << "match point at " << matchPoint << " for E = " << E << endl;
  	int matchPointIndex = int((matchPoint - xMin) / h);
  	matchPointIndex = max(4, matchPointIndex);
  	matchPointIndex = min(nPoints - 6, matchPointIndex);
  	// now we have to propagate left and right solutions until the matching point
  	solLR.clear();
    solLR.resize(matchPointIndex + 1);
    solLR[0].x = xMin;
    solLR[0].y = 0;
    solLR[1].x = xMin + h;
    solLR[1].y = 0.0001;  // this number should not be important, but it has to be small to avoid overflow

  	// cout << "matchPointIndex = " << matchPointIndex << " solLR.size() = " << solLR.size() << endl;
  	double h2 = h * h;
    for(int i = 2; i < solLR.size(); i++)
    {
        double x2  = xMin + (i - 2) * h;
        double x1  = xMin + (i - 1) * h;
        double x   = xMin + i * h;
        double p1  = 2 - 0.83333333333333 * h2 * shiftedPot(x1);
        double p2  = 1 + 0.08333333333333 * h2 * shiftedPot(x2);
        double p3  = 1 + 0.08333333333333 * h2 * shiftedPot(x);
        solLR[i].y = (p1 * solLR[i-1].y - p2 * solLR[i-2].y) / p3;
        solLR[i].x = x;
    }

  	solRL.clear();
  	solRL.resize(potential.size() - solLR.size());
    //fill the array, propagating the solRLution from the Right to the Left
    solRL[solRL.size() - 1].x = xMax;
    solRL[solRL.size() - 1].y = 0;
    solRL[solRL.size() - 2].x = xMax - h;
    solRL[solRL.size() - 2].y = 0.0001;

    for(int i = solRL.size() - 3; i > -1; i--)
    {
        int j = solRL.size() - 1 - i;  // distance to the end
        double x2  = xMax - (j - 2) * h;
        double x1  = xMax - (j - 1) * h;
        double x   = xMax - j * h;
        double p1  = 2 - 0.83333333333333 * h2 * shiftedPot(x1);
        double p2  = 1 + 0.08333333333333 * h2 * shiftedPot(x2);
        double p3  = 1 + 0.08333333333333 * h2 * shiftedPot(x);
        solRL[i].y = (p1 * solRL[i+1].y - p2 * solRL[i+2].y) / p3;
        solRL[i].x = x;
    }

  	// now we have to find the log derivative at the matching point
  	double v1 = solRL[0].y;
  	double d2 = ((25/12)*solLR[matchPointIndex].y-4*solLR[matchPointIndex-1].y+3*solLR[matchPointIndex-2].y-(4/3)* solLR[matchPointIndex-3].y+(1/4)* solLR[matchPointIndex-4].y);
  	double v2 = solLR[matchPointIndex].y;
  	double d1 = -((25/12)*solRL[0].y-4*solRL[1].y+3*solRL[2].y-(4/3)*solRL[3].y+(1/4) * solRL[4].y);

  	double logDer = (d1 * v2 - d2 * v1) / h;
  	//logDer /= (ld1 * ld1 + ld2 * ld2);
  	return logDer;
}

/**
 * Starting from the potential minima with a step of dEmin
 * computes diff(E), if a change of sign appears it means
 * that there is a zero thus store the range in the zeros vector.
 * These ranges are later used to find the root in a more refined phase.
 */
void Numerov::scanForZeroRanges(int nZeros) {
  // cout << "Scanning for zero ranges " << endl;
	findMinPot();
  // cout << "minimum of potential found" << endl;
	double E = minPot.y;
	zeros.clear();
	double lastDiff = diff(E), newDiff;
	while(zeros.size() < nZeros){
		newDiff = diff(E);
		if(newDiff * lastDiff < 0) {
			Range range(E - dEmin, E);
			zeros.push_back(range);
			 // cout << "zero in [" << range.eMin << ", " << range.eMax << "]" << endl;
		}
		lastDiff = newDiff;
		// just in case we hit a zero, which is not so unlikely to happen
		if(abs(lastDiff) < 1e-8) {
		    Range range(E - 1e-8, E + 1e-8);
		    zeros.push_back(range);
		    E += 1e-8;
		    //cout << "zero hit at " << E << endl;
		}

		E += dEmin;
	}
}

void Numerov::buildSol() {
	sol.clear();
	double scale = solRL[0].y / solLR[solLR.size() - 1].y;
	for (int i = 0; i < solLR.size(); i++) {
		Point p;
		p.x = solLR[i].x;
		p.y = scale * solLR[i].y;
		sol.push_back(p);
	}
	for (int i = 0; i < solRL.size(); i++) {
		sol.push_back(solRL[i]);
	}

	// normalize it
	double c = 0;
	for (int i = 1; i < sol.size(); i++)
		c += sol[i].y * sol[i].y * (sol[i].x - sol[i - 1].x);

	// TODO: this normalization process is not very good
	// users are suggested to re-normalize the outcome of this
	// by using some spline and quadrature adaptive algorithm for
	// numerical integration
	for (int i = 0; i < sol.size(); i++)
		sol[i].y /= sqrt(c);

	// we want to impose that the first maxima is alway positive
	// to avoid "jumps" while changing parameters in the potential
    auto derivative = [&](int i) {
        //https://en.wikipedia.org/wiki/Five-point_stencil
        return -sol[i + 2].y + 8 * sol[i + 1].y  - 8 * sol[i - 1].y + sol[i - 2].y;
    };

    double der = derivative(2);
    if(der < 0) // in this case change the overall sign
        for (int i = 0; i < sol.size(); i++)
            sol[i].y *= -1;
}

// given a solution this find its number
// of "nodes" which actually labels the solution
int Numerov::getSolIndex() {
    int mp = solLR.size();
    auto derivative = [&](int i) {
        //https://en.wikipedia.org/wiki/Five-point_stencil
        return -sol[i + 2].y + 8 * sol[i + 1].y  - 8 * sol[i - 1].y + sol[i - 2].y;
    };

    int flips = -1;
    //cout << "node found at ";
    double lastDer = derivative(2);
    for(int i = 3; i < sol.size() - 6; i++) {
        double newDer = derivative(i);

        if(lastDer * newDer < 0) {
            lastDer = newDer;
            // do not count the artifical nodes near the matching point
            if(mp - 50 < i && i < mp + 50)
                continue;
            //cout << sol[i].x << " ";
            flips++;
        }
    }
    //cout << endl;

    return flips;
}

void Numerov::findSpectrum(int nEigen){
  if(potential.size() < 2) {
    cout << "Please use the setPotential() function before using this one." << endl;
    return;
  }

  spectrum.clear();
  spectrum.potential = potential;
  scanForZeroRanges(nEigen);
  vector<Mode> modes;
  // check if a given index has been already computed
  auto hasBeenComputed = [&] (int index) {
    for (int i = 0; i < modes.size(); i++) {
      if(modes[i].index == index)
        return true;
    }

    return false;
  };

  const int MAX_DIVISIONS = 100;
  std::function<Range(double, double, int)> explore = [&](double from, double to, int N = 10) {
    //cout << "[DEBUG] exploring interval [" << from <<", " << to << ") N = " << N << endl;
    // this is too much, there should be something wrong (root not in the interval)
    if(N > MAX_DIVISIONS) {
      cout << "[WARN] The number of divisions " << N << " has gone beyond the limit for interval [" << from <<", " << to << "], aborting..." << endl;
      throw std::runtime_error("run time error");
    }

    // divide the interval in N and look for sign changes
    // if there is no success then attempts again with a finer grid
    double h = (to - from) / N;
    // find the first sign change
    Point lastVal(from, diffFunc(from));
    for(int i = 1; i < N; i++) {
      double p = from + i * h;
      Point newVal(p, diffFunc(p));
      if(lastVal.y * newVal.y < 0) {
        //cout << "New range for root [" << lastVal.x << ", "  << newVal.x << "]" << endl;
        Range r(lastVal.x, newVal.x);
        return r;
      }
      lastVal = newVal;
    }

    // if we reach this point we need to do a more refined search
    return explore(from, to, N + 50);
  };

  // find a solution in a range
  auto findSol = [&](Range r) {
    //cout << "finding sol in [" << r.eMin << ", " << r.eMax << "]..." << endl;
    double E = zbrent(diffFunc, r.eMin, r.eMax, tol);
    buildSol();
    int n = getSolIndex();
    //cout << "E(" << n << ") = " << E << endl;
    Mode mode(E, sol, n);
    return mode;
  };

  // using this for we get the right lowest nEigen eigenfunctions
  // most of the times, but sometimes we may skip some, see below.
  // int attempts = 0;
  // int const MAX_ATTEMPTS = 1;
  for (int i = 0; i < nEigen; i++) {
    Mode m = findSol(zeros[i]);
    modes.push_back(m);
    if(m.index >= nEigen) {
      //cout << "Possible bad index for eigenvalue " << m.index << endl;
      //cout << "Bound found, index " << m.index << endl;
      //break;
    }
  }

  auto findSpecSol = [&](int index) {
    int attempts = 0;
    while(!hasBeenComputed(index)) {
      int j = 0;
      // get the upper bound
      for(int i = 0; i < modes.size(); i++)
        if(modes[i].index > index) {
          j = i;
          break;
        }

        double lowLim = j == 0 ? minPot.y : modes[j - 1].energy;   // default is the lower limit, works for the ground state
        double minDiff = 1e-7;//(modes[i].energy - lowLim) / 100;
        // get the right interval to look for another root
        Range r = explore(lowLim + minDiff, modes[j].energy - minDiff, 20);
        Mode m = findSol(r);
        cout << "new mode found " << m.index << endl;
        modes.push_back(m);
        // we need to sort the modes and remove repeated after this insertion
        auto comp = [&](Mode m1, Mode m2) -> bool {
          return m1.energy < m2.energy;
        };
        auto rm = [&](Mode m1, Mode m2) -> bool {
          return m1.index == m2.index;
        };
        std::sort(modes.begin(), modes.end(), comp);
        std::unique(modes.begin(), modes.end(), rm);

        attempts++;
        if(attempts > 10) {
          cout << "[WARN] Too many attempts while finding eigenvalue " << index << ", computation is compromised." << endl;
          break;
        }
    }
  };

  // second recovery strategy
  // now we need to check if we get the right eigenfunctions
  // if we were asked to return the first 4 we don't want to return the 1,2,3 and 5!
  for (int i = 0; i < nEigen; i++)
    //TODO: finish this nice stuff
    if(false && !hasBeenComputed(i)) { // we need to look for a missing mode
      cout << "Looking back for E(" << i << ")..." << endl;
      try {
        findSpecSol(i);
      }catch(std::exception &e) {
        cout<<"Caught exception: "<<e.what()<<"\n";
        break;
      }
    }

    // finally safelly add all the modes found to the spectrum (already in a nice way)
    for (int i = 0; i < nEigen; i++)
      if(i < modes.size())
        spectrum.addMode(modes[i]);
      else{
        Mode m;
        spectrum.addMode(m);
      }
}
#endif
