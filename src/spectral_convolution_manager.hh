/**
 * @file   spectral_convolution_manager.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Nov 30 17:25:37 2015
 *
 * @brief  Class handling convolution of a field decomposed in several modes
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
#include <limits>
#include "convolution_manager.hh"
#include "ring_buffer.hh"
#include "data_register.hh"
/* -------------------------------------------------------------------------- */

class SpectralConvolutionManager : public ConvolutionManager {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  SpectralConvolutionManager(){};
  // Construct spectral manager with discretization used by the numerical integration,
  // a vector containing modal number (i.e. frequencies), the number of different fields and kernels
  SpectralConvolutionManager(Real dx, std::vector<Real> & modal_number,
			     UInt nb_fields, UInt nb_kernels);
  virtual ~SpectralConvolutionManager();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  
  // Initialize the manager given convolution cut and number of time steps
  // If blank is True, memory is not allocated (used to estimate the hardware capacity required to run a simulation)
  Idx init(Real cut, UInt nb_time_steps, bool blank);
  // (For a given kernel):
  // Use the default KernelFunctor and read pre-computed kernel values from file
  // Return poisson ratio read in the file for Real-check
  Real loadKernel(std::string filename, UInt kernel_id);
  // Get address of KernelFunctor of a given kernel
  // Should be used to set a different functor
  KernelFunctor ** getKernelFunctor(UInt kernel_id);
  // Store in buffer history an array of new values of size (2*nb_fields*nb_modes)
  // The scalar two stands for the real and imaginary parts of each value
  void storeFields(Real * new_vals);
  // Compute convolution between a given field and a given kernel
  // The result is an array of size (2*nb_modes)
  void computeConvolution(Real * res, UInt field_id,
			  UInt kernel_id);
  // Method used in restart framework of fields involved in convolutions (K and field_values)
  // pausing=true->generate restart files | pausing=false->restart simulation from existing files
  // If 3d simulation is restarted from 2d one, specify the number of modes (n_ele_fft={nele_x/2+1,nele_z})
  UInt restart(UInt side, bool pausing=false, std::vector<UInt> n_ele_fft={0,0});

private:

  virtual void initK();
  // Preintegrate nb_steps of convolution kernel (used mainly during restart) 
  void restartPreintegratedKernel(UInt nb_steps);
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // Number of modes handled by the spectral manager
  UInt nb_modes;
  // Number of fields
  UInt nb_fields;
  // Number of kernels
  UInt nb_kernels;
  // Number of values per field kept in memory
  Idx size_per_field;
  // Vector of modal numbers (i.e. frequencies)
  std::vector<Real> j_ksi;
  // Vector of convolution cut for each mode
  std::vector<UInt> mod_cut;
  // Position per mode of the memory share in K allocated for pre-integrated kernels 
  std::vector<UInt> k_start;
  // Vector of KernelFunctor used to compute the different kernels considered
  std::vector<KernelFunctor *> functors;
};

