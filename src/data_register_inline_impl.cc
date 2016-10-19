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
