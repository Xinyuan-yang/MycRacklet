/**
 * @file   simulation_driver_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Mon Oct 17 11:37:55 2016
 *
 * @brief  Implementation of the inline function of the SimulationDriver class
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
template<>
inline void SimulationDriver::setLoadingCase<_space_control>() {

  std::vector<UInt> n_ele = model.getNbElements();
  
  ctrl_loading.resize(3*(n_ele[0]-x_crack_start));
  for (UInt i=0; i < n_ele[0]-x_crack_start; ++i) {
    for (UInt j=0; j < 3; ++j) {
      ctrl_loading[3*i+j]= new_loading[j];
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void SimulationDriver::setLoadingCase<_time_control>() {

  Real ntim = model.getNbTimeSteps();
  
  if(ntim>0)
    ctrl_loading.reserve(3*ntim);
}

/* -------------------------------------------------------------------------- */
template<>
inline Real SimulationDriver::readSetLoadingCase<_time_control>(std::ifstream & file) {

  Real nb_time_steps;
  
  file >> nb_time_steps;
  
  ctrl_loading.resize(3*nb_time_steps);
  
  if(this->target_crack_speed)
    checkForTargetSpeed(file);
  
  Real line=0.0;
  
  for (UInt i =0; i < 3*nb_time_steps; ++i) { 
    if(file.good()){
    file>>line;   
    ctrl_loading[i]=line;
    }
    else
      cRacklet::error("!!! A problem occured with your loading file");
  }
  file.close();
  reading_time=-1;
  return nb_time_steps;

}

/* -------------------------------------------------------------------------- */
template<>
inline Real SimulationDriver::readSetLoadingCase<_space_control>(std::ifstream & file){

  Real nb_ele_file;
  std::vector<UInt> n_ele = model.getNbElements();
  
  file >> nb_ele_file;
  UInt file_jump;

  if (nb_ele_file<n_ele[0]) {
    std::stringstream err;
    err << "Unable to apply loading conditions !" << std::endl;
    err << "Your loading file was generated with only " << nb_ele_file <<" elements." << std::endl; 
    err << "Build your profile with at least " << n_ele[0] << " elements." << std::endl;
    cRacklet::error(err);
  } else {
    file_jump = nb_ele_file/n_ele[0];
  }
  ctrl_loading.resize(3*n_ele[0]);
  
  Real line=0.0;

  if(this->target_crack_speed)
    checkForTargetSpeed(file);
  
  for (UInt i =0; i < 3*n_ele[0]; ++i) { 
    if(file.good()) {
      for (UInt j = 0; j<file_jump; ++j){
	file>>line;   
      }
      ctrl_loading[i]=line;
    }
    else
      cRacklet::error("!!! A problem occured with your loading file");
  }
  file.close();
  return 0.;
}
