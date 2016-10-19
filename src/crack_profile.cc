/* -------------------------------------------------------------------------- */
#include "crack_profile.hh"
#include <vector>
#include <complex>
#include <fftw3.h>
#include <math.h>
/* -------------------------------------------------------------------------- */
CrackProfile CrackProfile::getStridedPart(UInt stride, UInt start) const {

  UInt total_lgth = heights.size()/stride;

  CrackProfile results(n,total_lgth);

  for (UInt i = 0; i < total_lgth; ++i) {
    results.heights[i] = heights[i*stride + start];
  }

 return results;

}

/* -------------------------------------------------------------------------- */
void CrackProfile::squareRoot() {

  for (UInt i = 0; i < heights.size(); ++i) {
    heights[i] = sqrt(heights[i]);
  }
}

/* -------------------------------------------------------------------------- */
void CrackProfile::stridedFFT(Real * output, UInt dim){
 
  int * N;
  fftw_complex * temp;
  UInt n_fft_per_dim = n[1]*(n[0]/2+1);

  temp = new fftw_complex[n_fft_per_dim*3]; 
  N = new int[2];

  N[0]=n[1];
  N[1]=n[0];

  fftw_plan p;

  p = fftw_plan_many_dft_r2c(2, N, dim, &heights[0], NULL, dim, 1, temp, NULL, 1, n_fft_per_dim, FFTW_ESTIMATE);
  
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  Real *it(output);

  for (UInt d = 0; d < dim; ++d) { 
    for (UInt i = 0; i < n_fft_per_dim-1; ++i) {
      for (UInt img = 0; img < 2; ++img) {
	*it = temp[d*n_fft_per_dim+i+1][img];
	++it;
      }
    }
  }
  delete[] temp;
  delete[] N;
}

/* -------------------------------------------------------------------------- */
void CrackProfile::backwardFFT(Real * input, UInt dim) {

  int * N;
  fftw_complex * temp;
  UInt n_fft_per_dim = n[1]*(n[0]/2+1);

  N = new int[2];

  N[0]=n[1];
  N[1]=n[0];
  
  temp = reinterpret_cast<fftw_complex*>(input);

  fftw_plan p;

  p = fftw_plan_many_dft_c2r(2, N, dim, temp, NULL, 1, n_fft_per_dim, &heights[0], NULL, dim, 1, FFTW_ESTIMATE);
  
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_cleanup();
  
  delete[] N;
}
