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
/* -------------------------------------------------------------------------- */


class CrackProfile {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  CrackProfile();
  CrackProfile(int size);
  virtual ~CrackProfile();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // resize the profile
  void SetGridSize(unsigned int a);
  // return a profile made of strided values
  CrackProfile getStridedPart(int stride, int start) const;
  // take the square root of each components
  void squareRoot();
  // compute a strided forward FFT and store only half spectrum without mean
  void stridedFFT(std::vector<double> & output, int stride);
  // compute a backward FFT from half spectrum (mean set as {0,0}) to a real strided profile
  void backwardFFT(std::complex<double> * input, int stride);
  // access to the i-th value in the array
  inline double & operator[](int i);
  inline const double & operator[](int i) const;
  // sum components by components
  inline CrackProfile operator+(const CrackProfile& q) const;
  // substraction components by components
  inline CrackProfile operator-(const CrackProfile& q) const;
  // multiplication components by components
  inline CrackProfile operator*(const CrackProfile& q) const;
  // multiplication each components with a real
  inline CrackProfile operator*(const double q) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // array of values for the profile itself
  std::vector<double> heights;
  // size of the array
  int n;
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
inline CrackProfile::CrackProfile(int a){
  SetGridSize(a);
}

/* -------------------------------------------------------------------------- */
inline void CrackProfile::SetGridSize(unsigned int a){
  heights.resize(a);
  n=a;
}

/* -------------------------------------------------------------------------- */
inline double & CrackProfile::operator[](int i){
  return heights[i];
};

/* -------------------------------------------------------------------------- */
inline const double & CrackProfile::operator[](int i) const{
  return heights[i];
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator+(const CrackProfile& q) const{
  
  CrackProfile results(n);

  for (int i = 0 ; i < n ; ++i){
    results.heights[i] = heights[i] + q.heights[i];
  }
  return results;
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator-(const CrackProfile& q) const{
  
  CrackProfile results(n);
  
  for (int i = 0 ; i < n ; ++i){
    results.heights[i] = heights[i] - q.heights[i];
  }
  return results;
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator*(const CrackProfile& q) const{
  
  CrackProfile results(n);
  
  for (int i = 0 ; i < n ; ++i){
    results.heights[i] = heights[i] * q.heights[i];
  }

  return results;
};

/* -------------------------------------------------------------------------- */
inline CrackProfile CrackProfile::operator*(const double q) const{
  
  CrackProfile results(n);
 
  for (int i = 0 ; i < n ; ++i){
    results.heights[i] = heights[i] * q;
  }
 
  return results;
};

#endif /* __CRACK_PROFILE_H__ */
