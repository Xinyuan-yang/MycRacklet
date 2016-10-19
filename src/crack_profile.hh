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

/* -------------------------------------------------------------------------- */
#ifndef __CRACK_PROFILE_H__
#define __CRACK_PROFILE_H__
/* -------------------------------------------------------------------------- */
#include <vector>
#include <fftw3.h>
#include <complex>
#include "cRacklet_common.hh"
/* -------------------------------------------------------------------------- */


class CrackProfile {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CrackProfile();
  CrackProfile(std::vector<UInt> size, UInt total_lgth);
  virtual ~CrackProfile();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // resize the profile
  void SetGridSize(std::vector<UInt> a, UInt dim);
  // get the total number of data stored
  UInt size() const {return heights.size();}
  // return a profile made of strided values
  CrackProfile getStridedPart(UInt stride, UInt start) const;
  // take the square root of each components
  void squareRoot();
  // compute a strided forward FFT and store only half spectrum without mean
  void stridedFFT(Real * output, UInt stride);
  // compute a backward FFT from half spectrum (mean set as {0,0}) to a real strided profile
  void backwardFFT(Real * input, UInt stride);
  // access to the i-th value in the array
  inline Real & operator[](UInt i);
  inline const Real & operator[](UInt i) const;
  // sum components by components
  inline CrackProfile operator+(const CrackProfile& q) const;
  // substraction components by components
  inline CrackProfile operator-(const CrackProfile& q) const;
  // multiplication components by components
  inline CrackProfile operator*(const CrackProfile& q) const;
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
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
inline CrackProfile::CrackProfile(){

}
/* -------------------------------------------------------------------------- */
inline CrackProfile::~CrackProfile(){

}
/* -------------------------------------------------------------------------- */
inline CrackProfile::CrackProfile(std::vector<UInt> a, UInt total_lgth){
  
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
inline CrackProfile CrackProfile::operator-(const CrackProfile& q) const{
  
  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] - q.heights[i];
  }
  return results;
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
inline CrackProfile CrackProfile::operator*(const Real q) const{

  UInt total_lgth = heights.size();

  CrackProfile results(n,total_lgth);
  
  for (UInt i = 0 ; i < total_lgth ; ++i){
    results.heights[i] = heights[i] * q;
  }
 
  return results;
};

#endif /* __CRACK_PROFILE_H__ */
