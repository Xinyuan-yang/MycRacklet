/**
 * @file   spectral_model.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Nov  9 19:10:51 2012
 *
 * @brief  Class performing the core operations of the spectral method
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

/* -------------------------------------------------------------------------- */
#ifndef __SPECTRAL_MODEL_H__
#define __SPECTRAL_MODEL_H__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "data_register.hh"
#include "crack_profile.hh"
#include "spectral_convolution_manager.hh"
#include "fracture_law.hh"
#include "contact_law.hh"
#include <vector>
#include <complex>
#include <math.h>
#include <iostream>
#include <fstream>
#ifdef CRACKLET_USE_LIBSURFER
#include "surface_generator_filter_fft.hh"
#include "surface_statistics.hh"
#endif
/* -------------------------------------------------------------------------- */

// Structure used to compute field variation across the interface
struct InterfaceFields{

  InterfaceFields(const std::vector<CrackProfile>* in_fields, UInt dimension);
  InterfaceFields(const CrackProfile* in_fields_top,
		  const CrackProfile* in_fields_bot, UInt dimension);
  // Initialization of fields
  void init(const CrackProfile* in_fields_top,
	    const CrackProfile* in_fields_bot, UInt dimension);
  // Compute delta_fields from fields values
  void computeJumpFields(); 
  // Interface fields[0] = top, fields[1] = bot
  std::vector<const CrackProfile*> fields;
  // Field jump through interface {spatial coordinates}
  CrackProfile fields_jump;
  // Field jump through interface {shear, normal}
  std::vector<CrackProfile> delta_fields;
  // Dimension
  UInt dim;

};

// Structure used to integrate energy values in time
struct Energetics{
  
  Energetics() {E=0; E_dot=0; E_dot_old=0;}
  // Compute E from E_dot and E_dot_old;
  void integrate(Real dt);

  // Energy dissipation
  Real E;
  // Energy dissipation rate (at t = it)
  Real E_dot;
  // Energy dissipation rate (at t = it-1)
  Real E_dot_old;
};

class SpectralModel : public DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  SpectralModel();
  SpectralModel(std::vector<UInt> nele, UInt nb_time_steps, std::vector<Real> dom_size, 
		Real nu_top, Real nu_bot, Real E_top, Real E_bot, 
		Real cs_top, Real cs_bot, UInt tcut_top, UInt tcut_bot, 
		bool overlap, UInt l_index, FractureLaw * fracturelaw, 
		ContactLaw * contactlaw, const std::string & simulation_summary,
		const std::string output_dir="./");

  virtual ~SpectralModel();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private :
  // Store pointer of the model fields into the data register 
  void registerModelFields();
  // Re-set beta parameters (the default value is set in cRacklet_common.hh)
  void resetStableTimeStep(Real new_beta);
  // Compute the modal numbers (i.e frequencies)
  template <UInt interface_dim>
  void initFrequency();
  // Initialize the two convolution managers (memory allocation, setting kernel functors)
  void initConvolutionManagers();
  // Compute the various convolution integrals
  template <UInt interface_dim>
  void computeConvolutions(Real * F_it, UInt side);
  // Compute normal velocities in case of relative slip
  void computeIndepNormalVelocities(UInt ix, UInt iz);
  // Compute shear velocities with a given shear strength
  void computeShearVelocities(Real strength, UInt elem);
  // Compute shear velocities in case of relative slip with a given strength
  void computeIndepShearVelocities(Real strength, UInt elem);
  // Compute velocities in the case of contact at crack type
  void computeContactVelocities(UInt ix, UInt iz);
 
public:

  // Initialization of the model. Specity a beta only to modify stable time step parameter
  // Default value defined in cRacklet_common.hh
  void initModel(Real beta=0.0);
  // Set a sinusoidal load distribution
  void sinusoidalLoading(Real min);
  // Set a brownian distributed loading
  // rms=root mean square, hurst=hurst exponent, q0=low cut_off, q1=roll_off, q2=high cut_off
  // !!! Required LibSurfer as an external library
  void brownianHeterogLoading(Real rms, long int seed, Real hurst, UInt q0,UInt q1, UInt q2);
  // Definition of the loading case
  void setLoadingCase(Real load_in, Real psi, Real phi);
  // Set loading case using a pre-computed loading file
  Real setLoadingCase(std::string loading_file, Real psi, Real phi);
  // update loading case
  void updateLoads();
  // Update point-wise loading conidtions using an uniform constant value per dimension  
  void updateLoads(Real * loading_per_dim);
  // update loading case from pre-computed loading condition
  UInt readUpdateLoads(Real start=0.0);
  // compute velocities at t=0
  void computeInitialVelocities();
  // Increases time step number
  void increaseTimeStep(){++it; displ_jump->computeJumpFields(); veloc_jump->computeJumpFields();};
  // Update displacements with velocities 
  void updateDisplacements();
  // Update material properties according to the related material laws
  void updateMaterialProp();
  // Compute FFT on displacement feelds
  void fftOnDisplacements();
  // Compute stresses convolution terms by Backward FFT
  void computeStress();
  // Compute velocities
  void computeVelocities();
  // Compute and energy release and energy release-rate at the current time step
  void computeEnergy();
  // dump the model paramters to simulation summary file
  void printSelfLoad(Real load, Real psi, Real phi);
  void printSelf();  

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  // Return the current simulation time
  Real getTime() {return it*beta*dxmin/X[0];}
  // Return the current simulation time step
  UInt getCurrentTimeStep() {return it;}
  // Return stable time step ratio beta
  Real getBeta() {return beta;}
  // Return plane discretization
  std::vector<Real> getElementSize() {return dx;}
  // Return the number of elements
  std::vector<UInt> getNbElements() {return n_ele;}
  // Return the number of time steps
  Real getNbTimeSteps() {return ntim;}
  // Return model dimension
  UInt getDim() {return dim;}
  // Return uniform loading vector used to set average interface loading conditions
  std::vector<Real> & getUniformLoading() {return uniform_loading;}
  // Return current crack tip position rightward from a starting position start
  UInt getCrackTipPosition(Real start);
  // Same but start expressed in terms of a discrete element number x_start
  UInt getCrackTipPosition(UInt x_start);
  // Return pointer to the FractureLaw
  FractureLaw * & getFractureLaw() {return fracture_law;}
  // Access to Energetic values
  const std::vector<Energetics> & getNormalDissipatedEnergy() {return E_nor;} 
  const std::vector<Energetics> & getShearDissipatedEnergy() {return E_shr;} 
  const std::vector<Energetics> & getFrictionalEnergy() {return E_fri;}
  const Energetics & getRadiatedEnergy() {return Eq;} 
  
public:
 
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */


private:
  
  // Number of elements must be a power of 2
  std::vector<UInt> n_ele;
  // Total number of elements (nex*nez)
  UInt total_n_ele;
  // Size of the domain
  std::vector<Real> X;
  //frequency 
  std::vector<Real> q00;
  //frequency k
  std::vector<Real> kk;
  //frequency m
  std::vector<Real> mm;
  // Number of time steps
  UInt ntim;
  // Poisson ratios of top and bottom material
  std::vector<Real> nu;
  // Top material wave speed
  Real cs_t;
  // Ratio of top and bottom shear wave spead
  Real ksi;
  // Top and bottom shear modulus
  std::vector<Real> mu;
  // Top and bottom density
  std::vector<Real> rho;
  // Time cut for tob and bottom domains
  std::vector<UInt> t_cut;
  // Overlapping tolerance (0 = no , 1 = yes)
  bool overlapping;
  // Index determining inside loading case for ind_crack>=index
  UInt load_index;
  // Associated fracture law
  FractureLaw * fracture_law;
  // Associated contact law
  ContactLaw * contact_law;
  // Dimension of interface
  UInt interface_dim;
  // Dimension of fields vectors (=3)
  UInt dim; 
  // Space step dx[0] = dx, dx[1] = dz (if 2d then dz is set to 0)
  std::vector<Real> dx;
  // minimum between dx and dz to compute dt
  Real dxmin;
  // Fundamental mode number
  std::vector<Real> q0;
  // Current time step
  UInt it;
  // number of spectral modes by direction for 3d nele_fft(nex, nez/2+1)
  std::vector<UInt> nele_fft;
  //total number of modes
  UInt total_nele_fft;
  // Displacements 
  std::vector<CrackProfile> displacements;  
  // Velocities
  std::vector<CrackProfile> velocities; 
  // Displacement jumps along the interface
  InterfaceFields * displ_jump;
  // Velocity jumps along the interface
  InterfaceFields * veloc_jump;
  // Stresses
  std::vector<CrackProfile> stresses;
  // Loads
  std::vector<CrackProfile> loads;
  // Point-wise ratio of loading compare to global average conditions
  Real * loading_ratio;
  // Tractions at the interface
  CrackProfile intfc_trac;
  // Normal Strength
  std::vector<Real> nor_strength;
  // Shear Strength
  std::vector<Real> shr_strength;
  // Cracking index just needed to better visualization of fracture process
  // Standard value: 0 = outside the crack, 1 = in the cohesive zone,
  // 2 = inside the crack, 3 = inside the contact zone
  // Other values can be defined arbitrarily at interface creation 
  std::vector<UInt> ind_crack;
  // Uniform loading condition per dimension applied along the interface
  std::vector<Real> uniform_loading;
  // Normal opening
  CrackProfile nor_opening;
  // Shear opening
  CrackProfile shr_opening;
  // Stable time step beta = cs*dt/dx
  Real beta;
  // Objects managing convolution integrals for top and bottom half-space
  SpectralConvolutionManager * convo_manager_top;
  SpectralConvolutionManager * convo_manager_bot;
  // Arrays allocated to store result of convolution integrals
  Real * h11u1,* h22u2,* h33u3,* h12u1,* h12u2,* h12u3,* h13u3,* h11u3,* h33u1;
  // Arrays allocated to store Fourier decomposition of the current displacements field
  Real * U_top, * U_bot;
  // Number of convolution kernel per domain used by the model
  UInt nb_kernels;
  // Ratio of top and bottom shear modulus
  Real zeta;
  // Dilatation over shear wave speed of top and bottom material
  std::vector<Real> eta;
  // Energy dissipated by normal opening
  std::vector<Energetics> E_nor;
  // Energy dissipated by shear opening
  std::vector<Energetics> E_shr;
  // Energy dissipated by friction 
  std::vector<Energetics> E_fri;
  // Radiated energy (Fukuyama 2005, Eq. 10)
  Energetics Eq;
  // Plan for the Fastest Fourier Transform in the West
  fftw_plan * plan;  
};

#include "spectral_model_inline_impl.cc"

#endif /* __SPECTRAL_MODEL_H__ */
