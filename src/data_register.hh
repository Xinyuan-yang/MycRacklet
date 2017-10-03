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
#ifndef __DATA_REGISTER_H__
#define __DATA_REGISTER_H__
/* -------------------------------------------------------------------------- */
#include "cRacklet_common.hh"
#include "crack_profile.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
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
  _normal_strength, //n_ele
  _shear_strength, //n_ele
  _frictional_strength, //n_ele
  _id_crack, //n_ele
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
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // List of grid points involved in the computations
  std::vector<UInt> index;
};

class DataRegister {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  DataRegister(){};
  virtual ~DataRegister(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  //Access to a given field returned in the type DataTypes
  static inline const DataTypes readData(DataFields my_field);
  //Register a computer object with a given name
  static void registerComputer(std::string computer_name, Computer * computer);
  //Access a registered computer object
  static Computer * getComputer(std::string computer_name);
  //Method returning current crack position searched between x_start and x_end (searching along z=0)
  static UInt getCrackTipPosition(UInt x_start, UInt x_end);
  
protected:
 
  //Initialize the resgister providing the output folder of simulation as well as its description
  void data_initialize(const std::string output_folder,const std::string description);
  //Finalize data register
  void data_finalize();
  //Register pointer to data associated with a given field name
  template<typename T>
  inline void registerData(DataFields my_field, T * in_data);
  //Launch all the registered computations for the given time
  void computeAll(Real time);
  
private:

  // Return the current daytime
  std::string getCurrentDaytime();

public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  // Ofstream of simulation summary file
  static std::ofstream out_summary;
  // Ofstream of simulation parameters file
  static std::ofstream out_parameters;
  // Output directory where simulation output will be written
  static std::string output_dir;

protected:
  // Map of registered data
  static std::map<DataFields,DataTypes> datas;
  // Map of registered computers
  static std::map<std::string,Computer*> computers;

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

#include "data_register_inline_impl.cc"

#endif /* __DATA_REGISTER_H__ */
