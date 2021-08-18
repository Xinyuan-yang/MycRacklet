/**
 * @file   rate_and_state_law.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Fri Feb 24 16:41:05 2017
 *
 * @brief  Class describing a rate and state friction law
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
#ifndef __RATE_AND_STATE_LAW__
#define __RATE_AND_STATE_LAW__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "interface_law.hh"
#include "crack_profile.hh"
#include <vector>
#include <math.h>
/* -------------------------------------------------------------------------- */
/**
 * @class  RateAndStateLaw rate_and_state_law.hh
 *
 * Class describing a rate and state friction law
 *
*/
class RateAndStateLaw : public InterfaceLaw {

  #include "rate_and_state_formulations.hh"
  
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:  
  /// Default constructor
  RateAndStateLaw()
    
    :InterfaceLaw() {

    UInt total_n_ele = n_ele[0]*n_ele[1];
   
    phi.resize(total_n_ele,0.0);
    cf.resize(total_n_ele,0.0);
    fric_strength.resize(total_n_ele);

    D.resize(total_n_ele);
    f_0.resize(total_n_ele);
    a.resize(total_n_ele);
    b.resize(total_n_ele);
    v_star.resize(total_n_ele);
    phi_star.resize(total_n_ele);
   
    this->registerData(_state_variable, &(phi));
    this->registerData(_friction_coefficient, &(cf));
    this->registerData(_frictional_strength, &fric_strength);
    this->registerData(_rands_D, &D);
    this->registerData(_rands_f_0, &f_0);
    this->registerData(_rands_a, &a);
    this->registerData(_rands_b, &b);
    this->registerData(_rands_v_star, &v_star);
    this->registerData(_rands_phi_star, &phi_star);
    
    sigma_0 = DataRegister::getParameter<Real>("sigma_0");

    shear_velo_jump = datas[_shear_velocity_jumps];
    dot_u_top = datas[_top_velocities];
    dot_u_bot = datas[_bottom_velocities];
    stresses = datas[_top_dynamic_stress];
    intfc_trac = datas[_interface_tractions];
    
    Real mu_top = this->getParameter<Real>("shear modulus top");
    Real mu_bot = this->getParameter<Real>("shear modulus bottom");
    Real cs_top = this->getParameter<Real>("shear wave speed top");
    Real cs_bot = this->getParameter<Real>("shear wave speed bottom");
    Real beta = this->getParameter<Real>("beta");
    Real dxmin = this->getParameter<Real>("delta min");
    
    if ((mu_top==mu_bot)&&(cs_top==cs_bot)) {
      c_s = cs_top;
      accoust = mu_top/c_s;
    }
    else
      cRacklet::error("Rate and state law is only implemented for homogeneous elastic properties");
    
    delta_t = beta*dxmin/c_s;
    V_0.resize(2*total_n_ele,0.);
  };

  /// Default destructor
  virtual ~RateAndStateLaw();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  /** Init the state evolution law
   */
  void initStateEvolution();
  /** Init Regularized state evolution
   */
  void initRegularizedStateEvolution(Real v0);
  /** Init the slip state evolution law
   */
  void initSlipStateEvolution();
  /** Init a given R&S formulation and state evolution
      standard rate and state formulation */
  void initStandardFormulation();
  /** Init formulation with a unique velocity-weakening branch
   */
  void initVelocityWeakeningFormulation();
  /** Init rate and state formulation with a regularized stick-to-slip transition
   */
  void initRegularizedFormulation(Real v0);
  /** Init rate and state formulation with a regularized stick-to-slip transition and weakening behavior at high velocities
   */
  void initRegularizedWeakeningFormulation(Real v0);
  /** Define the velocity prediction for each component (used before searching the initial steady state)
      @param v_0_pred : predictor of the velocity
   */ 
  void setVelocityPredictor(std::vector<Real> v_0_pred);
  /** Set the initial steady state velocity, if know a priori, for example in case of a restart. 
   */
  void setV0(std::vector<Real> v_0);
  /** Save the initial V0 profile
   */
  void saveV0();
  /** Initialize interface fields
   */
  void initInterfaceConditions();
  /** Compute the interface conditions but do not update them! (Used in pseudo-velocity driven systems)
   */
  std::vector<Real> computeNextAverageVelocity();
  /** Update the strength of the material in function of the opening profile
   */
  void updateInterfaceConditions();
  /** perturbe the state variable by adding \f$ \epsilon sin\left(2\pi k \frac{x}{X}\right) \f$
      @param epsilon : perturbation amplitude
      @param k : perturbation spatial frequency
   */
  void perturbState(Real epsilon, Real k);
  /** insert gaussian perturbation patch in the velocity field
   */
  void insertPerturbationPatch(std::vector<UInt> patch_limits, Real new_rate);
  /** insert gaussian perturbation in the velocity field 
   */
  void insertGaussianPerturbation(Real std_dev, Real amplitude);
  /** insert a skewed perturbation in the velocity field. The maximum of the perturbation is givenby the amplitude
   */
  void insertSkewedPerturbation(Real std_dev, Real amplitude, Real alpha);
  /** insert gaussian noise in the state field
   */
  void addGaussianNoiseToStateField(Real std_dev);
  /** insert perturbation from a given file to the state field
   */
  void insertPerturbationFromFile(std::string input_file);
  /** Method used in restart framework which still need to be implemented for the R&S law
   */
  void restart(bool pausing=false, UInt nele_2d=0);
    
private:
  /** initialize interface at an homogeneous the steady-state sliding conditions */ 
  void computeSteadyStateSliding();


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */

public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
   // Rate and state parameters
  std::vector<Real> D;
  std::vector<Real> f_0;
  std::vector<Real> a;
  std::vector<Real> b;
  std::vector<Real> v_star;
  std::vector<Real> phi_star;
  // State variable
  std::vector<Real> phi;
  // Friction coefficient
  std::vector<Real> cf;
  // Absract object respresenting the associated time evolution of the state variable
  std::shared_ptr<StateEvolution> state_evol;
  // Abstract object representing the associated rate and state formulation
  std::shared_ptr<RandSFormulation> formulation;
  // Frictional strengh
  std::vector<Real> fric_strength;
  // Vector of the initial steady-state sliding velocity
  std::vector<Real> V_0;
  
  CrackProfile * shear_velo_jump;
  CrackProfile * dot_u_top;
  CrackProfile * dot_u_bot;
  CrackProfile * stresses;
  CrackProfile * intfc_trac;
  
  Real sigma_0;
  Real c_s;
  Real accoust;
  Real delta_t;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

inline RateAndStateLaw::~RateAndStateLaw(){
  }

#endif /* __RATE_AND_STATE_LAW__ */
