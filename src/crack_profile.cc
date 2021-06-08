/**
 * @file   crack_profile.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Sun Jan  6 08:49:36 2013
 *
 * @brief  Implementation of the CrackProfile Class
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 * 
 * cRacklet is the result of a collaboration between the Computational Solid Mechanics 
 * Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland 
 * and the Department of Aerospace Engineering of the University of Illinois at 
 * Urbana-Champaign, United States of America.
 * 
 * cRacklet is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 * 
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with this program.  
 * If not, see <http://www.gnu.org/licenses/>.
 */
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
void CrackProfile::finalizeFFT() {
  
  if(data_fft) {
    fftw_destroy_plan(plan);
    delete[] data_fft;
  }
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
