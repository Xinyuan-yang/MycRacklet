/**
 * @file   simulation_driver.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Mon Oct 17 11:37:55 2016
 *
 * @brief  High level objects helping to drive SpectralModel
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef __SIMULATION_DRIVER__
#define __SIMULATION_DRIVER__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "spectral_model.hh"
/* -------------------------------------------------------------------------- */

// Enum to set the type of steady state crack speed algorithm
enum LoadControlType {  
  // To build load evolution as function of time
  _time_control,
  // To build load evolution as function of crack tip position
  _space_control
};

class SimulationDriver : public DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SimulationDriver(SpectralModel & model_to_drive): model(model_to_drive) {
    model.initModel();
  };
  SimulationDriver(SpectralModel & model_to_drive, Real target_speed,
		   Real crack_start);
  virtual ~SimulationDriver(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

private:
  // Change beta in case of fixed crack speed (integer number of time steps needed to break an element)
  Real adjustStableTimeStep();
  // Initialize vector needed to store loading evolution within steady state speed algorithm
  template<LoadControlType lc_type>
  void setLoadingCase();
  // Read loading evolution from ifstream to fix crack speed
  template <LoadControlType lc_type>
  Real readSetLoadingCase(std::ifstream & load_file);
  // Compute and store the evolution of loading conditions
  UInt runWritingStep();
  // Read and update the loading conditions
  UInt runReadingStep();
  // Solve one time step of SpectralModel in two phases
  void priorUpdates();
  void solveTimeStep();
  // Verify the compatibility of input file target speed with the one used to reset beta
  void checkForTargetSpeed(std::ifstream & file);
  // Algorithm tailoring loading condition for steady-state rupture. Return progression status (does crack_advance?)
  bool cstCrackFrontSpeed(UInt & x_tip);
public:

  // Init a simulation with constant loading conditions
  void initConstantLoading(Real cst_loading, Real psi, Real phi);
  // Init a simulation with evolving loading conditions following loading_file
  Real initLoadingFromFile(std::string loading_file, LoadControlType load_control=_time_control,
			   Real initial_loading=0., Real psi=0., Real phi=0.);  
  // Initialization before tailoring loading conditions to fix crack speed 
  void initConstantSpeed(Real initial_loading, Real psi, Real phi, Real average_crit_stress, 
			 Real spont_crack_length=0.0, LoadControlType load_control=_time_control);
  // Solve one time step of the simulation designed
  UInt solveStep();
  // Print tailored loading conditions to file
  void writeLoading(std::string load_file);
  // Launching artificially a through crack from position crack_start up to a given launched_size by artificially growing it at speed v_init*cs
  void launchCrack(Real crack_start, Real launched_size, Real v_init);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // Reference to the model to drive
  SpectralModel & model;
  // Type of algorithm used to control crack speed
  LoadControlType lc_type;
  //number of iteration needed to break one element at the desired crack speed 
  Real nb_t_by_elem;
  //Target crack speed (crack front speed/shear wave speed)
  Real target_crack_speed;
  //Center of the initial crack
  UInt x_crack_start;
  //Spontaneous crack length
  Real spont_crack;
  //Average critical stress used to normalized loading conditions
  Real av_crit_stress;
  //Loading imposed at the current time step
  std::vector<Real> new_loading;
  //Cohesive crack propagation  
  std::vector<Real> at;
  //Vector storing loading condition built to have a constant crack propagation speed 
  std::vector<Real> ctrl_loading;
  //Number of iteration needed to propagate to next grid spacing
  UInt countcfs;
  //Ratio: needed iteration over wanted iteration to propagate to next grid spacing
  std::vector<Real> atc;
  // Boolean to check when initialization phase is over
  bool c_s_initiation;
  // Time step at which the reading algorithm started (needed after the artifical growth of seed crack using launchCrack)
  int reading_time;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "simulation_driver_inline_impl.cc"

#endif /* __SIMULATION_DRIVER__ */
