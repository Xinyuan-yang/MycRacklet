/**
 * @file   data_register.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Aug 17 12:21:46 2016
 *
 * @brief  Mother class centralizing access to simulation data
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
#ifndef __DATA_REGISTER_H__
#define __DATA_REGISTER_H__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "crack_profile.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string.h>
/* -------------------------------------------------------------------------- */
//Identifier to help access to data wanted. The size of each field is indicated below
enum DataFields {
  _top_displacements, //n_ele*dim
  _bottom_displacements, //n_ele*dim
  _normal_displacement_jumps, //n_ele
  _shear_displacement_jumps, //n_ele
  _top_velocities, //n_ele*dim
  _bottom_velocities, //n_ele*dim
  _normal_velocity_jumps, //n_ele
  _shear_velocity_jumps, //n_ele
  _interface_tractions, //n_ele*dim
  _top_loading, //n_ele*dim
  _bottom_loading, //n_ele*dim
  _top_dynamic_stress, //n_ele*dim
  _bottom_dynamic_stress, //n_ele*dim
  
  _normal_strength, //n_ele
  _shear_strength, //n_ele
  _frictional_strength, //n_ele
  _id_crack, //n_ele
  _critical_normal_opening, //n_ele
  _critical_shear_opening, //n_ele

  _residual_normal_strength, //n_ele
  _residual_shear_strength, //n_ele
  
  // For viscoelastic law
  _lim_velocity, //n_ele
  
  _state_variable, //n_ele
  _friction_coefficient, //n_ele
  _rands_D, //n_ele
  _rands_f_0, //n_ele
  _rands_a, //n_ele
  _rands_b, //n_ele
  _rands_v_star, //n_ele
  _rands_phi_star, //n_ele
  
};

//Human readable names associated to each DataFields
static std::map<DataFields, std::string> datafields_name = {
  {_top_displacements, "Top displacements"},
  {_bottom_displacements, "Bottom displacements"},
  {_normal_displacement_jumps, "Normal displacement jumps"},
  {_shear_displacement_jumps, "Shear displacement jumps" },
  {_top_velocities, "Top velocities"},
  {_bottom_velocities, "Bottom velocities"},
  {_normal_velocity_jumps, "Normal velocity jumps"},
  {_shear_velocity_jumps, "Shear velocity jumps"},
  {_interface_tractions, "Interface tractions"},
  {_top_loading, "Top far-field loading"},
  {_bottom_loading, "Bottom far-field loading"},
  {_top_dynamic_stress, "Stress contribution from the history of the top displacements"},
  {_bottom_dynamic_stress, "Stress contribution from the history of the top displacements"},
  {_normal_strength, "Normal strength"},
  {_shear_strength, "Shear strength"},
  {_frictional_strength, "Frictional strength"},
  {_id_crack, "Crack index"},
  {_critical_normal_opening, "Critical normal opening displacement"},
  {_critical_shear_opening, "Critical shear opening displacement"},

  {_lim_velocity, "Limiting velocity for viscoelastic law"},

  {_state_variable, "State variable for R&S friction law"},
  {_friction_coefficient, "Coefficient of friction"},
  {_rands_D, "D coefficient of R&S friction law"},
  {_rands_f_0, "f_0 coefficient of R&S friction law"},
  {_rands_a, "a coefficient of R&S friction law"},
  {_rands_b, "b coefficient of R&S friction law"},
  {_rands_v_star, "v* coefficient of R&S friction law"},
  {_rands_phi_star, "phi* coefficient of R&S friction law"}
};

//Structure used to encapsulate different type of data pointers
//Only one of this pointer is non-NULL corresponding to the actual data type
struct DataTypes {
  //Initializing all pointer to NULL
  DataTypes(){
    vec_double=NULL;
    vec_uint=NULL;
    crack_prof=NULL;
  }
  //Operators converting a DataTypes object into its real pointer type
  template<typename T>
  inline operator T() const;
  template<typename T>
  inline operator T();
  //Converting a real pointer of type T into a generic DataTypes object
  template<typename T>
  inline void writeData(T*);

  //Members as list of currently available pointer type that can be converted into DataTypes  
  std::vector<Real> * vec_double;
  std::vector<UInt> * vec_uint;
  CrackProfile * crack_prof;
};

//Identifier to acces the different energy integrators implemented

enum IntegratorTypes {
  _shear_fracture_energy,
  _normal_fracture_energy,
  _frictional_energy,
  _radiated_energy,
};

// Virtual mother class for object performing computations on the model fields
class Computer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Computer(){};
  virtual ~Computer(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  // Virtual method launching computations at a given time step
  virtual void compute(Real time)=0;
  // Virtual method used in restart framework to save/retrieve computers values
  // pausing=true->save computers | pausing=false->retrieve computers values
  virtual void restart(std::fstream & file, bool pausing)=0;
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // List of grid points involved in the computations
  std::vector<UInt> index;
};

/**
 * @class DataRegister data_register.hh
 *
 * Mother class centralizing access to simulation data
 *
*/
class DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  /// Default Constructor
  DataRegister(){};
  /// Default Destructor
  virtual ~DataRegister(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  
  /// Access to a given field returned in the type DataTypes
  static inline const DataTypes readData(DataFields my_field);
  /** Register a global simulation parameter, whose type is templated T
      @param name (string) : name of the parameter
      @param value (T) : value
   */
  template<typename T>
  static void registerParameter(std::string name, T value);
  /** Register a set of global parameters from input file
      @param input_file (string):  path to the input file
   */
  static void readInputFile(std::string input_file);
  /** Check if a parameter is present before accessing it
      @param name (string) : name of the parameter
  */
  static inline bool hasParameter(std::string name);
  /** Get the value of a simulation parameter
      @param name (string) : name of the parameter
   */
  template<typename T>
  static inline T getParameter(std::string name);
  /// Register a computer object with a given name
  static void registerComputer(std::string computer_name, Computer * computer);
  /// Access a registered computer object
  static Computer * getComputer(std::string computer_name);
  /** Method returning current crack position searched between x_start and x_end (looking for the first point with state == 2)
      @param x_start (UInt) : index of the first element to start looking for the tip position
      @param x_end (UInt) : index of the last element to look for the tip position
      @param z_pos (UInt) : z index for 3D simulation. By default the position of the tip is investigated along z = 0
   */
  static UInt getCrackTipPosition(UInt x_start, UInt x_end, UInt z_pos=0);
  /// Method used in restart framework to load/export a vector from/to a data_file 
  /// pausing=true->generate restart data_file | pausing=false->load vector from existing data_file
  /// If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  template<typename T>
  static void restartData(std::vector<T> & my_data , const std::string data_file,
			  bool pausing, UInt nele_2d=0);  

  //Acessors
  
  /// Direct access to top velocities
  const CrackProfile * getTopVelocities();
  /// Direct access to bot velocities
  const CrackProfile * getBotVelocities();
  /// Direct access to shear velocity jumps
  const CrackProfile * getShearVelocityJumps();
  /// Direct access to normal velocity jumps
  const CrackProfile * getNormalVelocityJumps();
  /// Direct access to top displacement
  const CrackProfile * getTopDisplacements();
  /// Direct access to bot displacement
  const CrackProfile * getBotDisplacements();
  /// Direct access to shear displacement jumps
  const CrackProfile * getShearDisplacementJumps();
  /// Direct access to normal displacement jumps
  const CrackProfile * getNormalDisplacementJumps();
  /// Direct access to interface tractions
  const CrackProfile * getInterfaceTractions();

protected:
 
  //Initialize the register providing the output folder of simulation as well as its description
  void data_initialize(const std::string output_folder,const std::string description);
  //Finalize data register
  void data_finalize();
  //Register pointer to data associated with a given field name
  template<typename T>
  inline void registerData(DataFields my_field, T * in_data);
  //Launch all the registered computations for the given time
  void computeAll(Real time);
  // Method used in restart framework of fields saved in DataRegister
  // (displacements,velocities,crack_id,strength)
  // pausing=true->generate restart files | pausing=false->restart simulation from existing files
  // If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  static void restart(bool pausing=false, UInt nele_2d=0);
  
private:

  // Return the current daytime
  std::string getCurrentDaytime();
  // Method used in restart framework to save/retrieve computers values
  // pausing=true->save computers | pausing=false->retrieve computers values
  // If 3d simulation is restarted from 2d one, specify the number of elements along x (nele_2d=nele_x)
  static void restartComputer(bool pausing, UInt nele_2d);
  //Subrouting to restart 3d arrays from 2d
  template<typename T>
  static void restartDataFrom2d(std::vector<T> & my_data , std::fstream & file, UInt nele_2d);

  
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  /// Ofstream of simulation summary file
  static std::ofstream out_summary;
  /// Ofstream of simulation parameters file
  static std::ofstream out_parameters;
  /// Output directory where simulation outputs will be written
  static std::string output_dir;
  /// Output directory where restart outputs will be written
  static std::string restart_dir;

protected:

  // Map of registered data
  static std::map<DataFields,DataTypes> datas;
  // Map of registered computers
  static std::map<std::string,Computer*> computers;
  // Map to share global simulation variables
  static std::map<std::string,std::string> variables;
  // Dimension of fields vectors (=3)
  static UInt dim; 
  // Current time step
  static UInt it;
  // Number of elements must be a power of 2
  static std::vector<UInt> n_ele;
  // Stable time step beta = cs*dt/dx
  static Real beta;
  // Minimum between dx and dz used in the computation of compute dt
  static Real dxmin;
  // Ratio of top and bottom shear wave spead
  static Real ksi;
  // Dilatation over shear wave speed of top and bottom material
  static std::vector<Real> eta;
};

// Class computing integral over a defined sets of interface points (a.e. energy integration)
class Integrator : public Computer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  Integrator(){};
  // Construct an object by presicing the points, integrator type, starting time and area around point 
  Integrator(std::vector<UInt> integ_points, IntegratorTypes type, Real starting_time, Real dA)
  {this->index=integ_points; init(type,starting_time,dA);}
  virtual ~Integrator(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void init(IntegratorTypes type, Real starting_time, Real dA)
  {I=0; I_dot=0; I_dot_old=0; this->integ_type=type; t_old=starting_time; this->dA=dA;}
  
  // Integrate rate between last integration (t_old) and current time.
  void integrate(Real time);
  // Compute rate following the templated integrator types
  template<IntegratorTypes IT>
  inline void compute();
  
public:
  // Compute rate at a given time
  virtual inline void compute(Real time);
  // Method used in restart framework to save/retrieve integrators values
  // pausing=true->save computers | pausing=false->retrieve computers values
  virtual void restart(std::fstream & file, bool pausing) {
    if(pausing){
      file.write((char*) &(I), sizeof(Real));
      file.write((char*) &(I_dot_old), sizeof(Real));
      file.write((char*) &(t_old), sizeof(Real));
    }
    else {
      file.read((char*) &(I), sizeof(Real));
      file.read((char*) &(I_dot_old), sizeof(Real));
      file.read((char*) &(t_old), sizeof(Real));
    }
  }
  
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // Get integrated quantities
  Real getIntegration(){return I;}
  // Get current rate
  Real getRate(){return I_dot_old;}
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // Type of integration performed
  IntegratorTypes integ_type;
  // Integral value
  Real I;
  // Rate of change (at t = t_old+1)
  Real I_dot;
  // Rate of change (at t = t_old)
  Real I_dot_old;
  // Time of the last integration
  Real t_old;
  // Unit area of integration
  Real dA;
};

// Class computing integral over a given width following crack tip propagation.
// Crack propagation is tracked along z=0 between crack_start and crack_end
class SurfingIntegrator : public Integrator {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  SurfingIntegrator(UInt integration_width, IntegratorTypes type, Real starting_time, Real dA,
		    UInt crack_start, UInt crack_end) {
    init(type,starting_time,dA);
    this->integ_width = integration_width;
    this->crack_start = crack_start;
    this->crack_end = crack_end;
  }
  virtual ~SurfingIntegrator(){};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

private:
  // Update the list of integration points as function of current crack position
  void updateIntegrationPoints();
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // Compute rate at a given time
  inline void compute(Real time);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // Width of the integration domain surrounding crack tip positon
  UInt integ_width;
  // Crack tip position is tracked between crack_start and crack_end along z=0
  UInt crack_start;
  UInt crack_end;
  
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "data_register_inline_impl.hh"

#endif /* __DATA_REGISTER_H__ */
