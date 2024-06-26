/**
 * @file   data_register_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Aug 17 12:21:46 2016
 *
 * @brief  Implementation of the inline functions of DataRegister
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
inline const DataTypes DataRegister::readData(DataFields my_field) {
  return datas[my_field];
}

/* -------------------------------------------------------------------------- */
inline bool DataRegister::hasParameter(std::string name) {

  std::map<std::string,std::string>::iterator it = variables.find(name);

  if(it != variables.end())
    return true;
  else {
    return false;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T DataRegister::getParameter(std::string name) {

  std::map<std::string,std::string>::iterator it = variables.find(name);

  T val;

  if(it != variables.end()){
    std::stringstream stream(variables[name]);
    stream >> val;
    return val;
  }
  else {
    std::stringstream err;
    err << "No registered parameter named " << name << " !";
    cRacklet::error(err);
    return val;
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator const std::vector<UInt>*() const {
  return vec_uint;  
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator std::vector<UInt>*() {
  return vec_uint;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void DataTypes::writeData(std::vector<UInt> * in_data)  {
  vec_uint=in_data;  
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator const std::vector<Real>*() const {
  return vec_double;  
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator std::vector<Real>*() {
  return vec_double;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void DataTypes::writeData(std::vector<Real> * in_data)  {
  vec_double=in_data;  
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator const CrackProfile*() const {
  return crack_prof;  
}

/* -------------------------------------------------------------------------- */
template<>
inline DataTypes::operator CrackProfile*() {
  return crack_prof;  
}

/* -------------------------------------------------------------------------- */
template<>
inline void DataTypes::writeData(CrackProfile * in_data)  {
  crack_prof=in_data;  
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void DataRegister::registerData(DataFields my_field, T * in_data) {
  DataTypes new_data;
  new_data.writeData(in_data);
  datas[my_field] = new_data;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Integrator::compute<_shear_fracture_energy>() {

  const CrackProfile * shear_velo_jump = DataRegister::readData(_shear_velocity_jumps);
  const std::vector<Real> * shr_strength = DataRegister::readData(_shear_strength);

  for (UInt i = 0; i < this->index.size(); ++i) {
    this->I_dot += (*shr_strength)[index[i]] * fabs((*shear_velo_jump)[index[i]]);
    
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Integrator::compute<_normal_fracture_energy>() {

  const CrackProfile * normal_velo_jump = DataRegister::readData(_normal_velocity_jumps);
  const std::vector<Real> * nor_strength = DataRegister::readData(_normal_strength);

  for (UInt i = 0; i < this->index.size(); ++i) {
    this->I_dot += (*nor_strength)[index[i]] * fabs((*normal_velo_jump)[index[i]]);
    
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Integrator::compute<_frictional_energy>() {

  const CrackProfile * shear_velo_jump = DataRegister::readData(_shear_velocity_jumps);
  const std::vector<Real> * fric_strength = DataRegister::readData(_frictional_strength);
  const std::vector<UInt> * id_crack = DataRegister::readData(_id_crack);
  for (UInt i = 0; i < this->index.size(); ++i) {
    if ((*id_crack)[index[i]]==2){
      this->I_dot += (*fric_strength)[index[i]] * fabs((*shear_velo_jump)[index[i]]);}
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Integrator::compute<_radiated_energy>() {

  const CrackProfile * interface_tractions = DataRegister::readData(_interface_tractions);
  const CrackProfile * sigma_0_top = DataRegister::readData(_top_loading);
  const CrackProfile * sigma_0_bot = DataRegister::readData(_bottom_loading);
  const CrackProfile * v_top = DataRegister::readData(_top_velocities);
  const CrackProfile * v_bot = DataRegister::readData(_bottom_velocities);
  
  for (UInt i = 0; i < this->index.size(); ++i) {
    for (UInt d = 0; d < 3; ++d) {
      UInt idx = 3*index[i]+d;
      this->I_dot += ((*sigma_0_top)[idx]-(*interface_tractions)[idx])*((*v_top)[idx])
	- ((*sigma_0_bot)[idx]-(*interface_tractions)[idx])*((*v_bot)[idx]);
    }
  }  
}

/* -------------------------------------------------------------------------- */
inline void Integrator::compute(Real time) {

  switch (integ_type) {
    
  case _shear_fracture_energy:
    compute<_shear_fracture_energy>();
    break;
  case _normal_fracture_energy:
    compute<_normal_fracture_energy>();
    break;
  case _frictional_energy:
    compute<_frictional_energy>();
    break;
  case _radiated_energy:
    compute<_radiated_energy>();
    break;
    
  default:
    std::stringstream err;
    err << "*** IntegratorType (" << integ_type 
	<< ") is not implemented in data_register_inline_impl.hh" << std::endl;
    cRacklet::error(err);
    break;
  }

  this->I_dot *= dA;
  integrate(time);
}

/* -------------------------------------------------------------------------- */
inline void SurfingIntegrator::compute(Real time) {

  updateIntegrationPoints();
  Integrator::compute(time);
}
