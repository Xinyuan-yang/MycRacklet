/**
 * @file   spectral_model.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Nov  9 19:10:51 2012
 *
 * @brief  Class processing spectral method for the study of crack propagation at interfaces
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

#include "crack_profile.hh"
#include "circular_buffer.hh"
#include "fracture_law.hh"
#include "contact_law.hh"
#include <vector>
#include <fftw3.h>
#include <complex>
#include <math.h>

struct InterfaceFields{

  InterfaceFields(const std::vector<CrackProfile> * in_fields, int dimension)
  { fields = in_fields; dim = dimension;}
  // Compute delta_fileds from fields values
  void computeJumpFields();
 
  // Interface fields
  const std::vector<CrackProfile> * fields;
  // Field jump through interface {shear, normal}
  std::vector<CrackProfile> delta_fields;
  // Dimension
  int dim;

};

struct Energetics{
  
  Energetics() {E=0; E_dot=0; E_dot_old=0;}
  // Compute E from E_dot and E_dot_old;
  void integrate(double dt);

  // Energy dissipation
  double E;
  // Energy dissipation rate (at t = it)
  double E_dot;
  // Energy dissipation rate (at t = it-1)
  double E_dot_old;
};

class SpectralModel {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  SpectralModel();
  SpectralModel(int nele, int nb_time_steps, double dom_size, double crack_size, 
		double nu_top, double nu_bot, double E_top, double E_bot, 
		double cs_top, double cs_bot, int tcut_top, int tcut_bot, 
		bool overlap, unsigned int l_index, FractureLaw * fracturelaw, 
		ContactLaw * contactlaw);

  virtual ~SpectralModel();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
private :
  // Reader of a kernel binary files generated by inverse_serial.f
  void readKernel(std::string f90file, int ik);
  // Templated binary reader
  template<class Bed>
  void readBinary(std::ifstream & file, Bed & val);
  // Fonction that returns kernel H33 = J1(x)/x, with J1 = Bessel function
  // of the 1st kind
  double H33(double x);
  // interpolate kernel values from pre-computed discrete loaded kernel 
  double interpolateFromKernel(double val, int index);
  // compute normal velocities in case of relative slip
  void computeIndepNormalVelocities(int elem);
  // compute shear velocities with a given shear strength
  void computeShearVelocities(double strength, int elem);
  // compute shear velocities in case of relative slip with a given strength
  void computeIndepShearVelocities(double strength, int elem);
  // compute velocities in the case of contact at crack type
  void computeContactVelocities(int elem);
 
public:

  // Initialization of the model
  void initModel();
  // Definition of the loading case
  void setLoadingCase(double load_in, double load_out, double psi, double phi);
   // update loading case
  void updateLoads();
  // compute velocities at t=0
  void computeInitialVelocities();
  // Increases time step number
  void increaseTimeStep(){ ++it;};
  // Update displacements with velocities 
  void updateDisplacements();
  // Update material properties according to the related material laws
  void updateMaterialProp();
  // Increase by beta an additionnal loading of amplitude ampl defined to trigger slip at position trigger
  void updateTriggerLoading(double alpha, double beta, double ampl);
  // Load pre-computed convolution kernel from files
  void computeKernels();
  // Pre-integrate convolution kernels
  void preintegratedKernels();
  // Compute FFT on displacement feelds
  void fftOnDisplacements();
  // Compute convolution integrals
  void computeTimeConvolution();
  // Compute stresses convolution terms by Backward FFT
  void computeStresses();
  // Compute velocities
  void computeVelocities();
  // Compute and print energy released and energy release-rate at the current time step
  void computeEnergy();
  // Generate a pulse (x delta) and length lx at the center of the domain
  void generatePulse(double delta, double lx);
  // dump the model paramters in a given ofstream
  void printSelf(std::ofstream & parameters_file, std::ofstream & summary);
   
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  // Reference to displacements profiles
  const std::vector<CrackProfile> & getDisplacements() {return displacements;}
  // Reference to velocities profiles
  const std::vector<CrackProfile> & getVelocities() {return velocities;} 
  // Reference to interface tractions profile
  const CrackProfile & getInterfaceTractions() {return intfc_trac;}
  // Reference to strength arrays
  std::vector<double> & getNormalStrength() {return nor_strength;} 
  std::vector<double> & getShearStrength() {return shr_strength;} 
  // Reference to cracking index array
  std::vector<unsigned int> & getCrackingIndex() {return ind_crack;}
  // Access to Energetic values
  const std::vector<Energetics> & getNormalDissipatedEnergy() {return E_nor;} 
  const std::vector<Energetics> & getShearDissipatedEnergy() {return E_shr;} 
  const std::vector<Energetics> & getFrictionalEnergy() {return E_fri;} 
  // Return the main simualtion parameters of the model
  const std::vector<double> getSimulationsParameters() { 
    return {X, (double)n_ele, beta, (double)dim, a_0};}

public:
 
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  
  // Number of elements must be a power of 2
  int n_ele;
  // Size of the domain
  double X;
  // Initial crack size
  double a_0;
  // Number of time steps
  int ntim;
  // Poisson ratios of top and bottom material
  std::vector<double> nu;
  // Top material wave speed
  double cs_t;
  // Ratio of top and bottom shear wave spead
  double ksi;
  // Top and bottom shear modulus
  std::vector<double> mu;
  // Time cut for tob and bottom domains
  std::vector<int> t_cut;
  // Overlapping tolerance (0 = no , 1 = yes)
  bool overlapping;
  // Index determining inside loading case for ind_crack>=index
  unsigned int load_index;
  // Associated fracture law
  FractureLaw * fracture_law;
  // Associated contact law
  ContactLaw * contact_law;
  // Dimension of displacements vector (=3)
  int dim; 
  // Space step
  double dx;
  // Fundamental mode number
  double q0;
  // Current time step
  int it;
  // Size of the half domain
  int nele_fft;
  // Displacements 
  std::vector<CrackProfile> displacements;  
  // Velocities
  std::vector<CrackProfile> velocities; 
  // Objects to compute interface jumps
  InterfaceFields * displ_jump;
  InterfaceFields * veloc_jump;
  // Stresses
  std::vector<CrackProfile> stresses;
  // Loads
  std::vector<CrackProfile> loads;
  // Tractions at the interface
  CrackProfile intfc_trac;
  // Normal Strength
  std::vector<double> nor_strength;
  // Shear Strength
  std::vector<double> shr_strength;
  // Cracking index (0 = outside the crack, 1 = in the cohesive zone,
  // 2 = inside the crack, 3 = inside the contact zone, 
  // 4 = high fracture toughness region
  std::vector<unsigned int> ind_crack;
  // Loads applied inside the crack zone
  std::vector<double> loads_insd;
  // Loads applied outside the crack zone
  std::vector<double> loads_outsd;
  // Normal opening
  CrackProfile nor_opening;
  // Shear opening
  CrackProfile shr_opening;
  // Stable time step beta = cs*dt/dx
  double beta;
  // Pre-computed convolution kernel values
  std::vector<std::vector<double> > kernel;
  // Number of convolution kernel per domain used by the model
  int nb_kernels;
  // Time cut for each mode for top and bottom domains
  std::vector<std::vector<int> > t_cut_j;
  // Sampling spacing used for the different kernel precomputations
  std::vector<double> delta_smpl;
  // Pré-integrated convolution kernel
  std::vector<std::vector<std::vector<double> > > K;
  // Displacements Fourier coefficients
  std::vector<CircularBuffer<double> > U;
  // Integrated convolutions convo = [H11U1, H12U2, H12U1, H22U2, H33U3]
  std::vector<std::vector<std::complex<double> > > convo;
  // Ratio of top and bottom shear modulus
  double zeta;
  // Dilatation over shear wave speed of top and bottom material
  std::vector<double> eta;
  // Energy dissipated by normal opening
  std::vector<Energetics> E_nor;
  // Energy dissipated by shear opening
  std::vector<Energetics> E_shr;
  // Energy dissipated by friction 
  std::vector<Energetics> E_fri;
  // Plan for the Fastest Fourier Transform in the West
  fftw_plan * plan;

};

#include "spectral_model_impl.cc"

#endif /* __SPECTRAL_MODEL_H__ */
