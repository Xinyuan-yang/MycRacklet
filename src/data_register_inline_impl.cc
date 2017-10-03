inline const DataTypes DataRegister::readData(DataFields my_field) {
  return datas[my_field];
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
    if ((*id_crack)[index[i]]==3)
      this->I_dot += (*fric_strength)[index[i]] * fabs((*shear_velo_jump)[index[i]]);
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
	<< ") is not implemented in data_register_inline_impl.cc" << std::endl;
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
