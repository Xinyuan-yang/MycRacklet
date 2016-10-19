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
  static const DataTypes readData(DataFields my_field);
  
protected:
 
  //Initialize the resgister providing the output folder of simulation as well as its description
  void data_initialize(const std::string output_folder,const std::string description);
  //Finalize data register
  void data_finalize();
  //Register pointer to data associated with a given field name
  template<typename T>
  inline void registerData(DataFields my_field, T * in_data); 
  
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
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "data_register_inline_impl.cc"

#endif /* __DATA_REGISTER_H__ */
