#ifndef _CHEBSPEC_
#define _CHEBSPEC_
#include <stdlib.h>
#include "solvspec.hpp"
// arpack++ can be installed very easily as described in
// https://github.com/m-reuter/arpackpp/blob/master/INSTALL.md
// #include "areig.h"
// #include "ardnsmat.h"

extern "C" {
  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  // product C= alphaA.B + betaC
  void dgemm_(char* TRANSA, char* TRANSB, const int* M,
                const int* N, const int* K, double* alpha, double* A,
                const int* LDA, double* B, const int* LDB, double* beta,
                double* C, const int* LDC);
  // diagonalizer non symmetric real matrices
  void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
              int* lda, double* wr, double* wi, double* vl, int* ldvl,
              double* vr, int* ldvr, double* work, int* lwork, int* info );

}

class ChebSpec : public SolvSpec {
	private:
	  int N = -1;
	  int L = 0;
	  double b, a, scal;
	  vector<double> x;
	  vector<double> V;
	  vector<vector<double>> Tm;
	  vector<vector<double>> TmInt;
	  vector<vector<double>> D;
	  vector<vector<double>> D2;
	  vector<vector<double>> Ehat;
	  vector<vector<double>> UE;
	  vector<vector<double>> US;
	  void invMat(double*, int);
	  void showMatrix(double* A, int Nr);
	  virtual void setN(int n);

	public:
		virtual void findSpectrum(int nEigen);
		virtual void setPotential(vector<Point>);
	  ChebSpec() {
	    cout << "Initializing ChebSpec" << endl;
	    setN(200);
	  }
};

void ChebSpec::setN(int n) {
  if(N == n)
    return;
  N = n;
  L = n + 1;
  cout << "computing chebyshev matrices, N = " << N << endl;
  // initialize x
  x.resize(L);
  V.resize(L);
  for(int i = 0; i < L; i++)
    x[i] = -cos(PI * i / N);

  // Tm, the matrix of the Chebyshev polynomial values in the grid. Checked.
  Tm.resize(L);
  for(int i = 0; i < L; i++) {
    Tm[i].resize(L);
    for(int j = 0; j < L; j++) {
      if(i == 0)      // first
        Tm[i][j] = 1;
      else if(i == 1) // second
        Tm[i][j] = x[j];
      else
        Tm[i][j] = 2 * x[j] * Tm[i - 1][j] - Tm[i - 2][j];
    }
  }

  // TmIntInt, the matrix of the integrated from -1 Chebyshev polynomial values in the grid. Checked.
  TmInt.resize(L);
  for(int i = 0; i < L; i++) {
    TmInt[i].resize(L);
    for(int j = 0; j < L; j++) {
      if(i == 0)      // first
        TmInt[i][j] = x[j] + 1;
      else if(i == 1) // second
        TmInt[i][j] =  0.5 * (x[j] * x[j] - 1);
      else if(i < L - 1)
        TmInt[i][j] = Tm[i + 1][j] / (2 * (i + 1)) - Tm[i - 1][j] / (2 * (i - 1)) + pow(-1, i + 1) / (pow(i, 2) - 1);
      else   // case i = L - 1, right extreme
        TmInt[i][j] = (2 * x[j] * Tm[i][j] - Tm[i - 1][j]) / (2 * (i + 1)) - Tm[i - 1][j] / (2 * (i - 1)) + pow(-1, i + 1) / (pow(i, 2) - 1);
    }
  }

  // D matrix: the integral from -1 to x. Checked.
  D.resize(L);
  for(int i = 0; i < L; i++) {
    D[i].resize(L);
    for(int j = 0; j < L; j++) {
      double s = 0;
      for(int k = 0; k < L; k++) {
        double f = 1;
        if(k == 0 || k == L - 1) // prime '' sum
          f = 0.5;

        s += f * Tm[k][j] * TmInt[k][i];
      }

      D[i][j] = (2. / N) * s;
      // prime '' sum
      if(j == 0 || j == L - 1)
        D[i][j] = 0.5 * D[i][j];

      // cout << D[i][j] << " ";
    }
    // cout << endl;
  }

  // finally the D2 = D*D matrix
  // check, the sum of the last row has to be 2 and 0 the sum of the first
  D2.resize(L);
  for(int i = 0; i < L; i++) {
    D2[i].resize(L);
    for(int j = 0; j < L; j++) {
      double s = 0;
      for(int k = 0; k < L; k++)
        s += D[i][k] * D[k][j];

      D2[i][j] = s;
    }
  }

  // check, the sum of the last row has to be 2
  double s = 0;
  for(int k = 0; k < L; k++)
    s += D2[0][k];

  cout << "sum D2_Nk " << s;
}

void ChebSpec::setPotential(vector<Point> p) {
  SolvSpec::setPotential(p);
  a = potential[0].x;
  b = potential[potential.size() - 1].x;
  scal = pow((b - a) / 2, 2);
  cout << "Interval [" << a << ", " << b << "] mapped to [-1,1]" << endl;
}

void ChebSpec::findSpectrum(int nEigen) {
  if(potential.size() < 2) {
    cout << "Please use the setPotential() function before using this one." << endl;
    return;
  }

  // compute V
  for(int i = 0; i < L; i++)
    V[i] = scal * potFunc(0.5 * ((b - a) * x[i] + b + a));


  // compute Ehat
  Ehat.resize(L);
  for(int i = 0; i < L; i++) {
    Ehat[i].resize(L);
    for(int j = 0; j < L; j++) {
      Ehat[i][j] = D2[i][j] * V[j];
    }
  }

  cout << "[DEBUG] compute UE and US" << endl;
  // compute UE and US
  UE.resize(L);
  US.resize(L);
  for(int i = 0; i < L; i++) {
    UE[i].resize(L);
    US[i].resize(L);
    for(int j = 0; j < L; j++) {
      UE[i][j] = 0.5 * (x[i] + 1) * D2[N][j];
      US[i][j] = 0.5 * (x[i] + 1) * Ehat[N][j];
    }
  }

  cout << endl << "[DEBUG] compute A and B " << endl;
  // compute the A and B matrices, remove the first/last row/column
  int Nr = N - 1;
  double A[Nr * Nr], B[Nr * Nr], B1A[Nr * Nr];

  //lapack routines use the column-major order!
  for(int i = 0; i < Nr; i++) {
    for(int j = 0; j < Nr; j++) {
      // define A
      A[i + j * Nr] = -US[i + 1][j + 1] + Ehat[i + 1][j + 1];
      if(i == j)   // add the identity matrix
        A[i + j * Nr] = A[i + j * Nr] - 1;
      // now B
      B[i + j * Nr] = D2[i + 1][j + 1] - UE[i + 1][j + 1];
    }
  }
  // we need to find the eigenvalues of A x = lambda B x
  invMat(B, Nr);      //B is now its inverse
  //showMatrix(B, Nr);
  char TRANS = 'N';
  double ALPHA = 1.0;
  double BETA = 0.0;
  double wr[Nr], w[Nr], wi[Nr], vl[Nr*Nr], vr[Nr*Nr];
  //Calculate B^-1 * A = B1A
  dgemm_(&TRANS, &TRANS, &Nr, &Nr, &Nr, &ALPHA, B, &Nr, A, &Nr, &BETA, B1A, &Nr);
  // diagonalize B1A
  cout << "[DEBUG] diagonalize B1A...";
  int n = Nr, lda = Nr, ldvl = Nr, ldvr = Nr, info, lwork;
  double wkopt;
  double* work;
  /* Query and allocate the optimal workspace */
  lwork = -1;
  char vec[] = "Vectors";
  dgeev_( vec, vec, &n, B1A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Solve eigenproblem */
  dgeev_( "Vectors", "Vectors", &n, B1A, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
         work, &lwork, &info );
  cout << "done" << endl;

  // DEBUG: print eigenvalues
  for(int i = 0; i < Nr; i++)
    cout << "E[" << i << "] = " << wr[i] / scal << endl;
  // ENDDEBUG

  free( (void*) work );
  //delete[] WORK;
}

void ChebSpec::showMatrix(double* A, int Nr) {
  cout << endl;
  for(int i = 0; i < Nr; i++) {
    for(int j = 0; j < Nr; j++) {
      cout << A[i + j * Nr] << " ";
    }
    cout << endl;
  }
}

void ChebSpec::invMat(double* A, int Nr)
{
    int *IPIV = new int[Nr+1];
    int LWORK = Nr*Nr;
    double *WORK = new double[LWORK];
    int INFO;
    dgetrf_(&Nr,&Nr,A,&Nr,IPIV,&INFO);
    dgetri_(&Nr,A,&Nr,IPIV,WORK,&LWORK,&INFO);
    delete[] IPIV;
    delete[] WORK;
}

#endif
