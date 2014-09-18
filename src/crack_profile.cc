/* -------------------------------------------------------------------------- */
#include "crack_profile.hh"
#include <vector>
#include <fftw3.h>
#include <math.h>
#include <complex>
/* -------------------------------------------------------------------------- */
CrackProfile CrackProfile::getStridedPart(int stride, int start) const {

  int size = n/stride;

  CrackProfile results(size);

for (int i = 0; i < size; ++i) {
  results.heights[i] = heights[i*stride + start];
 }

 return results;

}

/* -------------------------------------------------------------------------- */
void CrackProfile::squareRoot() {

  for (int i = 0; i < n; ++i) {
    heights[i] = sqrt(heights[i]);
  }
}

/* -------------------------------------------------------------------------- */
void CrackProfile::stridedFFT(std::vector<double> & output, int dim){
 
  int nele_fft = n/2; 
  int * N;
  fftw_complex * temp;
  int n_per_dim = n/dim;
  int n_fft_per_dim = n_per_dim/2+1;

  temp = new fftw_complex[nele_fft+3];

  N = new int[dim];

  for (int i = 0; i < dim; ++i) {
    N[i]=n_per_dim;
  }
 
  fftw_plan p;

  p = fftw_plan_many_dft_r2c(1, N, dim, &heights[0], NULL, dim, 1, temp, NULL, 1, n_fft_per_dim, FFTW_ESTIMATE);
  
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  for (int d = 0; d < dim; ++d) { 
    for (int i = 0; i < n_fft_per_dim-1; ++i) {
      for (int img = 0; img < 2; ++img) {

	output[2*d*(n_fft_per_dim-1)+2*i+img] = temp[d*n_fft_per_dim+i+1][img];
    
      }
    }
  }
  delete[] temp;
  delete[] N;
}

/* -------------------------------------------------------------------------- */
void CrackProfile::backwardFFT(std::complex<double> * input, int dim) {

  int nele_fft = n/2;
  int * N;
  fftw_complex * temp;
  int n_per_dim = n/dim;
  int n_fft_per_dim = n_per_dim/2+1;

  N = new int[dim];

  for (int i = 0; i < dim; ++i) {
    N[i]=n_per_dim;
  }

  for (int i = 0; i < nele_fft; ++i) {
  
    temp = reinterpret_cast<fftw_complex*>(input);

  }

  fftw_plan p;

  p = fftw_plan_many_dft_c2r(1, N, dim, temp, NULL, 1, n_fft_per_dim, &heights[0], NULL, dim, 1, FFTW_ESTIMATE);
  
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  delete[] N;
}

