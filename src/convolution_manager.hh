/**
 * @file   convolution_manager.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Nov 13 12:15:34 2015
 *
 * @brief  Class handling convolution integral between given kernel(s) and field(s)
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
#ifndef __CONVOLUTION_MANAGER__
#define __CONVOLUTION_MANAGER__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "ring_buffer.hh"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h> 
#include <gsl/gsl_sf_bessel.h>
/* -------------------------------------------------------------------------- */
class ConvolutionManager {

public:
  // Functor use to compute kernel 
  struct KernelFunctor{
    KernelFunctor(){};
    virtual ~KernelFunctor(){};
    // Compute kernel value as function of tau
    // Default behavior interpolate the result from pre-computed values 
    virtual Real operator()(Real tau);
    // Load pre-computed values from a file
    Real loadKernelFromFile(std::string filename);
  private:
    // Templated binary reader use to read pre-computed kernel files
    template<class Bed>
    void readBinary(std::ifstream & file, Bed & val);
  private:
    // Vector containing pre-computed values
    std::vector<Real> kernel;
    // Pre-computing sampling spacing used when generating kernel files
    Real delta_smpl; 
  };
  // Functor computing kernel H33 = J1(x)/x, 
  // with J1 = Bessel function of the 1st kind
  struct H_33 : public KernelFunctor{
    H_33(){};
    virtual ~H_33(){};
    Real operator()(Real tau){

      if(tau!=0)
	return gsl_sf_bessel_J1(tau)/tau;
      else
	return 0.5;
    }
  };
   
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  
  ConvolutionManager(){};
  // Construct manager with a given discretization used by the numerical integration
  ConvolutionManager(Real dx){this->dx=dx;this->func=NULL;this->K=NULL;this->field_values=NULL;}
  virtual ~ConvolutionManager();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // Initialize the manager with a given convolution cut
  void init(UInt cut);
  // Set the KernelFunctor used to compute kernel values
  void loadKernel(KernelFunctor * func);
  // Use the default KernelFunctor and read pre-computed kernel values from file
  // Return poisson ratio read in the file for double-check
  Real loadKernel(std::string filename);
  // Store a new field value in buffer history 
  inline void storeFields(Real new_val);
  // Compute convolution between field and pre-integrated convolution kernel
  inline Real computeConvolution();
  // Print nb_points of kernel function between from and to
  void dumpKernel(std::string filename, Real from, Real to, UInt nb_points);

protected:
  // Allocate the array memory used to store field history
  void allocateMemory(Idx size);
  // Initialize a given RingBuffer 
  // using an array starting at a given address with a given size of allocated memory  
  void initRingBuffer(RingBuffer<Real> & buffer_to_init, 
		      Real * buffer_start, UInt buffer_size);
  // Initialize and allocate memory for the array containing pre-integrated kernel
  virtual void initK();
  // Store a field value in a given RingBuffer
  inline void storeFields(RingBuffer<Real> & buffer_destination,
			  Real new_val);
  // Pre-integrate convolution kernel at time i using a given KernelFunctor
  // Pre-integration is multiplied by a scalar alpha and stored in memory share starting at K+k_start
  inline void preintegrateKernel(UInt i, KernelFunctor * funct,
				 UInt k_start=0, Real alpha=1.);
  // Compute convolution between a pre-integrated kernel stored in memory share starting at K+k_start
  // a field managed by a given buffer with a given convolution cut
  inline Real convolute(RingBuffer<Real> & buffer, 
			  UInt cut, UInt k_start=0);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  // Array where pre-integrated kernel values are stored
  Real * K;
  // Array where history of field value is stored
  Real * field_values;
 // Ring-buffer object(s) manipulating field_values with FIFO operations
  RingBuffer<Real> * field;
  // Disretization used for numerical integration
  Real dx;
  // Functor used to compute the convolution kernel
  KernelFunctor * func;
  // Number of field values kept in memory   
  Idx size;
  
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "convolution_manager_inline_impl.hh"

#endif /* __CONVOLUTION_MANAGER__ */
