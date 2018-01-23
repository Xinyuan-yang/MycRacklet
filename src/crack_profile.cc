/* -------------------------------------------------------------------------- */
#include "crack_profile.hh"
#include <vector>
#include <complex>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#if defined (_OPENMP)
#include <omp.h>
#endif
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
void CrackProfile::initFFT(bool forward_fft, UInt dim) {

  if(heights.size()==0)
    cRacklet::error("CrackProfile with size 0 cannot be initialized");
  
  int * N;
  N = new int[2];
  
  N[0]=n[1];
  N[1]=n[0];

  UInt n_fft_per_dim = n[1]*(n[0]/2+1);
  
  data_fft = new fftw_complex[n_fft_per_dim*3];
 
  if (forward_fft)
    plan = fftw_plan_many_dft_r2c(2, N, dim, &heights[0], NULL, dim, 1, data_fft, NULL, 1, n_fft_per_dim, FFTW_ESTIMATE);
  else
    plan = fftw_plan_many_dft_c2r(2, N, dim, data_fft, NULL, 1, n_fft_per_dim, &heights[0], NULL, dim, 1, FFTW_ESTIMATE);

  delete[] N;
  
}

/* -------------------------------------------------------------------------- */
void CrackProfile::squareRoot() {

  for (UInt i = 0; i < heights.size(); ++i) {
    heights[i] = sqrt(heights[i]);
  }
}

/* -------------------------------------------------------------------------- */
void CrackProfile::stridedFFT(Real * output, UInt dim){
  
  fftw_execute(this->plan);
   
  Real *it(output);
  UInt n_fft_per_dim = n[1]*(n[0]/2+1);
  
  for (UInt d = 0; d < dim; ++d) { 
    for (UInt i = 0; i < n_fft_per_dim-1; ++i) {
      for (UInt img = 0; img < 2; ++img) {
	*it = data_fft[d*n_fft_per_dim+i+1][img];
	++it;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void CrackProfile::backwardFFT(Real * input, UInt dim) {

  UInt n_fft_per_dim = n[1]*(n[0]/2+1);

  memcpy(data_fft, input, sizeof(Real)*3*n_fft_per_dim*2);

  fftw_execute(this->plan);
}
