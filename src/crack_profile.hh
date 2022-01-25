/**
 * @file   crack_profile.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @date   Sun Jan  6 08:49:36 2013
 *
 * @brief  Class dealing with real paramters along a crack-profile
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
#ifndef __CRACK_PROFILE_H__
#define __CRACK_PROFILE_H__
/* -------------------------------------------------------------------------- */
#include <vector>
#include <fftw3.h>
#include <complex>
#include <algorithm>
#include "cRacklet_common.hh"
#if defined (_OPENMP)
#include <omp.h>
#endif
/* -------------------------------------------------------------------------- */


class CrackProfile {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CrackProfile();
  CrackProfile(std::vector<UInt> size, UInt total_lgth);
  virtual ~CrackProfile(){};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // resize the profile
  void SetGridSize(std::vector<UInt> a, UInt dim);
  // initialize fftw library objects (only needed if fft transforms are wanted)
  // crack_profile in initialized to perform eihter forward (forward_fft=true)
  // or backward ((forward_fft=false) operations
  void initFFT(bool forward_fft, UInt dim);
  // finalize and free fftw library objects. Note that cleanup is only call once by the DataRegister
  void finalizeFFT();
  // get the total number of data stored
  UInt size() const {return heights.size();}
  // get the shape of stored data
  std::vector<UInt> shape() const {return n;}
  // return a profile made of strided values
  CrackProfile getStridedPart(UInt stride, UInt start) const;
  // return the value of maximum height
  Real getMaxValue() const {return *(max_element(heights.begin(), heights.end()));}
  // take the square root of each components
  void squareRoot();
  // compute a strided forward FFT and store only half spectrum without mean
  void stridedFFT(Real * output, UInt stride);
  // compute a backward FFT from half spectrum (mean set as {0,0}) to a real strided profile
  void backwardFFT(Real * input, UInt stride);
  // return the array of values
  std::vector<Real> & getValues(){return heights;}
  const std::vector<Real> & getValues() const {return heights;}
  // access to the i-th value in the array
  inline Real & operator[](UInt i);
  inline const Real & operator[](UInt i) const;
  // sum components by components
  inline CrackProfile operator+(const CrackProfile& q) const;
  inline void operator+=(const CrackProfile& q);
  // substraction components by components
  inline CrackProfile operator-(const CrackProfile& q) const;
  inline void operator-=(const CrackProfile& q);
  // multiplication components by components
  inline CrackProfile operator*(const CrackProfile& q) const;
  inline void operator*=(const Real q);
  // multiplication each components with a real
  inline CrackProfile operator*(const Real q) const;
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // array of values for the profile itself
  std::vector<Real> heights;
  // number of elements in each grid direction
  std::vector<UInt> n; //{nex,nez}
  // fftw library complex data (output of forward fft or input of backward fft)
  fftw_complex * data_fft;
  // fftw library object handling the transformations
  fftw_plan plan;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline CrackProfile::CrackProfile(){

  data_fft=NULL;
}

/* -------------------------------------------------------------------------- */
inline CrackProfile::CrackProfile(std::vector<UInt> a, UInt total_lgth) : CrackProfile() {
  
  n = a;
  heights.resize(total_lgth);
}

/* -------------------------------------------------------------------------- */
inline void CrackProfile::SetGridSize(std::vector<UInt> a, UInt dim){
  
  UInt total_n_ele=1;

  for (UInt i = 0; i < a.size(); ++i) {
    total_n_ele *= a[i]; 
  }

  n = a;

  heights.resize(total_n_ele*dim);
  
}

/* -------------------------------------------------------------------------- */
inline Real & CrackProfile::operator[](UInt i){
  return heights[i];
};

/* -------------------------------------------------------------------------- */
inline const Real & CrackProfile::operator[](UInt i) const{
  return heights[i];
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator+(const CrackProfile& q) const{
  
  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] + q.heights[i];
  }
  return results;
};

/* -------------------------------------------------------------------------- */
inline void CrackProfile::operator+=(const CrackProfile& q) {

  UInt total_size = heights.size();
  
  if(q.heights.size()!=total_size)
    cRacklet::error("CrackProfile should be of same size for += operation");
  else{  
    for (UInt i = 0 ; i < total_size; ++i){
      heights[i] = heights[i] + q.heights[i];
    }
  }
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator-(const CrackProfile& q) const{
  
  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] - q.heights[i];
  }
  return results;
};

/* -------------------------------------------------------------------------- */
inline void CrackProfile::operator-=(const CrackProfile& q) {

  UInt total_size = heights.size();
  
  if(q.heights.size()!=total_size)
    cRacklet::error("CrackProfile should be of same size for -= operation");
  else{  
    for (UInt i = 0 ; i < total_size; ++i){
      heights[i] = heights[i] - q.heights[i];
    }
  }
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator*(const CrackProfile& q) const{
   
  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] * q.heights[i];
  }

  return results;
};

/* -------------------------------------------------------------------------- */
inline void CrackProfile::operator*=(const Real q) {

  UInt total_size = heights.size();
  
  for (UInt i = 0 ; i < total_size; ++i){
    heights[i] = heights[i]*q;
  }
};


/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator*(const Real q) const{

  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] * q;
  }
 
  return results;
};

#endif /* __CRACK_PROFILE_H__ */
