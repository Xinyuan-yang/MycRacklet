/**
 * @file   interfacer.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Sep 11 13:17:04 2014
 *
 * @brief This class is used to set interface properties  
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __INTERFACER__
#define __INTERFACER__
/* -------------------------------------------------------------------------- */
#include "data_register.hh"
#include "spectral_model.hh"
#include "fracture_law.hh"
#ifdef CRACKLET_USE_LIBSURFER
#include "surface_generator_filter_fft.hh"
#include "surface_statistics.hh"
#endif
/* -------------------------------------------------------------------------- */
enum FractureLawType {
  _linear_coupled_cohesive
};

template<FractureLawType>
class Interfacer : public DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Interfacer(SpectralModel & model) : fracture_law(model.getFractureLaw()) {

    shr_strength = datas[_shear_strength];
    nor_strength = datas[_normal_strength];
    ind_crack = datas[_id_crack];
    
    dx.resize(2);
    dx[0] = model.getElementSize()[0];
    dx[1] = model.getElementSize()[1];
    n_ele.resize(2);
    n_ele[0] = model.getNbElements()[0];
    n_ele[1] = model.getNbElements()[1];
    total_n_ele = n_ele[0]*n_ele[1];
  }

  virtual ~Interfacer(){};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  //Some lexical conventions:
  //THROUGH: designate an z-invariant area
  //POLAR ASPERITY: designate an asperity made of weaker then stronger interface properties
  // polarity = 0: strong-weak / polarity = 1: weak-strong

  void applyInterfaceCreation();
  // create a uniform layer on the entire interface
  void createUniformInterface(Real crit_nor_opening, Real crit_shr_opening, 
			      Real max_nor_strength, Real max_shr_strength);
  // Tune interface properties from a text file starting at a given position
  void insertPatternfromFile(std::string filename, UInt origin=0);
  // create an heterogeneous interface following normal distribution of strength
  void createNormalDistributedInterface(Real crit_nor_opening, 
					Real crit_shr_opening, 
					Real max_nor_strength, 
					Real max_shr_strength,
					Real stddev, Real seed);
  // create an heterogeneous interface with a brownian distribution of strength
  // rms=root mean square, hurst=hurst exponent, q0=low cut_off, q1=roll_off, q2=high cut_off
  // !!! Required LibSurfer as an external library
  void createBrownianHeterogInterface(Real crit_nor_opening, 
				      Real crit_shr_opening, 
				      Real max_nor_strength, 
				      Real max_shr_strength,
				      Real rms, long int seed,
				      Real hurst=0.8, UInt q0=4,
				      UInt q1=4, UInt q2=32);

  // create a z-invariant(="through") area between x=start and x=end of given cracking_index and
  // with properties given by
  // new_prop = ratio*current_prop, if variation_rather_than_ratio=0
  // new_prop = ratio+current_prop, if variation_rather_than_ratio=1
  void createThroughArea(Real area_start, Real area_end,
			 UInt cracking_index,
			 Real ratio_max_nor_strength=1., 
			 Real ratio_max_shr_strength=1.,
			 Real ratio_crit_nor_opening=1., 
			 Real ratio_crit_shr_opening=1.,
			 bool variation_rather_than_ratio=0);
  // create a crack between x=crack_start and x=crack_end) 
  void createThroughCrack(Real crack_start, Real crack_end);
  // create a wall between x=wall_start and x=wall_end
  void createThroughWall(Real wall_start, Real wall_end);
  // create an asperity made of weaker(-delta) then stronger(+delta) areas at given position and given width
  void createThroughPolarAsperity(Real position, Real width,
				  Real delta_max_nor_strength, 
				  Real delta_max_shr_strength,
				  Real delta_crit_nor_opening, 
				  Real delta_crit_shr_opening, 
				  bool polarity);
  // create an interface made of multiple weaker(-delta) and stronger(+delta) areas
  // A given number of asperities is inserted between x=start and x=end
  // Return effective x_end position (accounting that Real end is rounded to match discretization)
  UInt createThroughMultiPolAsperity(Real start, Real end,
				     Real number,
				     Real delta_max_nor_strength, 
				     Real delta_max_shr_strength,
				     Real delta_crit_nor_opening, 
				     Real delta_crit_shr_opening, 
				     bool polarity);

  // create an interface with a centered crack.
  void createThroughCenteredCrack(Real initial_crack_size, Real crit_nor_opening, Real crit_shr_opening, 
				  Real max_nor_strength, Real max_shr_strength);
  // create an interface with a left-sided crack.
  void createThroughLeftSidedCrack(Real initial_crack_size, Real crit_nor_opening, Real crit_shr_opening, 
				  Real max_nor_strength, Real max_shr_strength);
 
  // create right propagating through crack meeting a circular asperity whose
  // strength = ratio_strength*interface_strength and crit_opening = ratio_crit_open*interface_crit_opening
  void createRightPropagatingCrackRoundAsp(Real initial_crack_size, Real crit_nor_opening,
					   Real crit_shr_opening, Real max_nor_strength,
					   Real max_shr_strength, Real radius,
					   std::vector<Real> asp_ctr, Real ratio_strength, Real ratio_crit_open=1.);
  
  // create an interface without initial cohesion between top and bottom solids
  void createIncohIntfc();
 
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // Fracture law object to build
  FractureLaw ** fracture_law;
  // Shear and normal strength arrays of the interface
  std::vector<Real> * shr_strength;
  std::vector<Real> * nor_strength;
  std::vector<Real> crit_n_open;
  std::vector<Real> crit_s_open;
  // Cracking index
  std::vector<UInt> * ind_crack;
  // Number of elements at the interface
  std::vector<UInt> n_ele;
  // Total number of elements
  UInt total_n_ele;
  // Space step
  std::vector<Real> dx; 
  // Interface dimension 1 if 2D space and 2 if 3D space
  UInt interface_dim;

};

#include "interfacer_inline_impl.cc"

#endif  /* __INTERFACER__ */
