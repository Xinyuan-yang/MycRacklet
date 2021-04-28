/**
 * @file   simulation_driver.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
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

/**
 * @class  SimulationDriver simulation_driver.hh
 *
 * High level objects helping to drive SpectralModel
 *
*/
class SimulationDriver : public DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /** Constructor for simulation driver
      @param model_to_drive : (&SpectralModel) Reference to the model to drive
      @param beta : (Real) Pre-factor of the stable time step. Default value in cRacklet_commom.hh
  */
  SimulationDriver(SpectralModel & model_to_drive, Real beta=0.0): model(model_to_drive) {
    model.initModel(beta);
  };
  /** Constructor for simulation driver with target driving speed
      @param model_to_drive : (&SpectralModel) Reference to the model to drive
      @param target_speed : (Real) Target speed of the rupture, in fraction of the shear wave speed (<1)
      @param crack_start : (Real) Starting position of the crack to drive
  */
  SimulationDriver(SpectralModel & model_to_drive, Real target_speed,
		   Real crack_start);
  /// Standard destructor
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
  void solveTimeStep();
  // Verify the compatibility of input file target speed with the one used to reset beta
  void checkForTargetSpeed(std::ifstream & file);
  // Algorithm tailoring loading condition for steady-state rupture. Return progression status (does crack_advance?)
  bool cstCrackFrontSpeed(UInt & x_tip);

public:
  /** Init a simulation with constant loading conditions
      @param cst_loading : (Real) Norm of the loading vector
      @param psi : (Real) Angle of the loading with respect to the x axis
      @param phi : (Real) Angle of the loading with respect to the y axis
  */
  void initConstantLoading(Real cst_loading, Real psi, Real phi);
  /** Init a simulation with evolving loading conditions following loading_file
      @param loading_file : (string) File containing the loading
      @param load_control : (LoadControlType) Either _time_control (loading is linked to time) or _space_control (loading is associated to a crack position
      @param initial_loading : (Real) Initial guess for the loading, default 0
      @param psi : (Real) Angle of the loading with respect to the x axis, default 0
      @param phi : (Real) Angle of the loading with respect to the y axis, default 0
  */
  Real initLoadingFromFile(std::string loading_file, LoadControlType load_control=_time_control,
			   Real initial_loading=0., Real psi=0., Real phi=0.);  
  /** Initialization before tailoring loading conditions to fix crack speed
      @param initial_loading : (Real) Initial value of the loading vector
      @param psi : (Real) Angle of the loading with respect to the x axis
      @param phi : (Real) Angle of the loading with respect to the y axis
      @param average_max_stress : (Real) Maximum possible value for the remote loading
      @param spont_crack_length : (Real) spontaneous crack length
      @param load_control : (LoadControlType) Either _time_control (loading is linked to time) or _space_control (loading is associated to a crack position
      @param max_average_stress : (Real) Upper bound of the applied far-field loading. Can be defined as a fraction of average_max_stress
      @param griffith_length : (Real) largest static crack size used to define the minimum applied far-field loading
  */
  void initConstantSpeed(Real initial_loading, Real psi, Real phi, Real average_max_stress, 
			 Real spont_crack_length=0.0, LoadControlType load_control=_time_control,
			 Real load_upper_bound=0.9, Real griffith_length=0.);
  /// Solve one time step of the simulation
  UInt solveStep();
  /** Print tailored loading conditions to file
      @param load_file: (string) name of the load file
   */
  void writeLoading(std::string load_file);
  /** Launching artificially a through crack from position crack_start up to a given launched_size by artificially growing it at speed v_init*cs
      @param crack_start: (Real) initial crack position 
      @param lauched_size: (Real) end position of the target crack      
      @param v_init: (string) growth velocity, multiplied by cs the shear wave speed
      @param one_side_propagation: (bool) one side propagation, by default is True
   */
  void launchCrack(Real crack_start, Real launched_size, Real v_init, bool one_side_propagation=true);

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
  //Average maximal stress from the cohesive law used to normalized loading conditions
  Real av_max_stress;
  //Upper and lower bounds of the far-field loading (used in case of steady-state simulation
  Real max_load;
  Real min_load;
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

#include "simulation_driver_inline_impl.hh"

#endif /* __SIMULATION_DRIVER__ */
