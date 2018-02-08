// DEPRECATED: for some reason it randomly stops in my computer
// taken from https://personal.broadinstitute.org/anshul/projects/chromatinVariation/src/ase_cpp/deps/vcftools/dgeev.h
#include <cmath>
#include <stdlib.h>

extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
		       int *lda, double *wr, double *wi, double *vl,
		       int *ldvl, double *vr, int *ldvr,
		       double *work, int *lwork, int *info);

void dgeev(double *a, int n, double *Er, double *Ei, double **Evecs);
void dgeev_ftoc(double *in, double **out, int rows, int cols);
void dgeev_sort(double *Er, double *Ei, double **Evecs, int N);

void dgeev(double *a, int n, double *Er, double *Ei, double **Evecs)
{
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  double *vl, *vr, *work;

  jobvl = 'N';
  jobvr = 'V';
  lda = n;

  ldvl = n;
  vl = new double[n*n];
  ldvr = n;
  vr = new double[n*n];
  lwork = -1;
  double wkopt;
  dgeev_(&jobvl, &jobvr, &n, a, &lda, Er, Ei, vl, &ldvl, vr, &ldvr, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  dgeev_(&jobvl, &jobvr, &n, a, &lda, Er, Ei, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
  dgeev_ftoc(vr, Evecs, n, ldvr);
  dgeev_sort(Er, Ei, Evecs, n);

  delete [] a;
  delete [] vl;
  delete [] vr;
  free( (void*) work );
}


double* dgeev_ctof(double **in, int rows, int cols)
{
  double *out;
  int i, j;

  out = new double[rows*cols];
  for (i=0; i<rows; i++)
	  for (j=0; j<cols; j++)
		  out[i+j*cols] = in[i][j];
  return(out);
}


void dgeev_ftoc(double *in, double **out, int rows, int cols)
{
	int i, j;

	for (i=0; i<rows; i++)
		for (j=0; j<cols; j++)
			out[i][j] = in[i+j*cols];
}

void dgeev_sort(double *Er, double *Ei, double **Evecs, int N)
{
	double temp, *E2;
	int i, j, k;

	E2 = new double[N];
	for (i=0; i<N; i++)
		E2[i] = Er[i]*Er[i]+Ei[i]*Ei[i];

	for (j=0; j<N; j++)
		for (i=0; i<N-1; i++)
			if (fabs(E2[i])>fabs(E2[i+1]))
			{
				temp = E2[i]; E2[i] = E2[i+1]; E2[i+1] = temp;
				temp = Er[i]; Er[i] = Er[i+1]; Er[i+1] = temp;
				temp = Ei[i]; Ei[i] = Ei[i+1]; Ei[i+1] = temp;

				for (k=0; k<N; k++)
				{
					temp = Evecs[k][i];
					Evecs[k][i] = Evecs[k][i+1];
					Evecs[k][i+1] = temp;
				}
			}

	delete [] E2;
}
