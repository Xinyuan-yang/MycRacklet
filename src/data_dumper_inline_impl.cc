inline void Dumper::dump() {
  
  switch (field_name) {

  case _top_displacements:
  case _bottom_displacements:
  case _normal_displacement_jumps:
  case _shear_displacement_jumps:
  case _top_velocities:
  case _bottom_velocities:
  case _normal_velocity_jumps:
  case _shear_velocity_jumps:
  case _interface_tractions:
  case _top_dynamic_stress:
  case _bottom_dynamic_stress:
  case _top_loading:
  case _bottom_loading:
    dump<CrackProfile>();
    break;

  case _normal_strength:
  case _shear_strength:
  case _frictional_strength:
  case _state_variable:
  case _friction_coefficient:
    dump<std::vector<Real> >();
    break;
    
  case _id_crack:
    dump<std::vector<UInt> >();
    break;
  default:
    std::stringstream err;
    err << "*** The DataTypes to dump (" << field_name 
	<< ") is not defined in data_dumper.cc" << std::endl;
    cRacklet::error(err);
    break;
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline void Dumper::dump_binary(Bed & data, UInt size) {

  std::streamsize bed_size = sizeof(data[0]);
  
  for (UInt i = 0; i < size; ++i) {
    file.write((char*)(&(data[i*stride+start])),bed_size);
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline void Dumper::dump_text(Bed & data, UInt size) {

  for (UInt i = 0; i < size; ++i) {
    file << std::scientific << std::setprecision(9) << data[i*stride+start] << " ";
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline void PointsDumper::dump_binary(Bed & data, UInt size) {

  for (UInt i = 0; i < size; ++i) {
    Real to_write = (Real)(data[start*size+i]);
    file.write((char*)(&to_write),sizeof(Real));
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline void PointsDumper::dump_text(Bed & data, UInt size) {

  for (UInt i = 0; i < size; ++i) {
    file << std::scientific << std::setprecision(9) << data[start*size+i] << " ";
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline void Dumper::dump() {

  const Bed * data_ptr = DataRegister::readData(field_name);
  const Bed & data = *(data_ptr);
  UInt size = (UInt)(ratio*data.size());
 
  switch (output_format) {

  case _text:
    if(PointsDumper * ptr = dynamic_cast<PointsDumper*>(this))
      ptr->dump_text(data,size);
    else
      this->dump_text(data,size);
    break;
  case _binary:
    if(PointsDumper * ptr = dynamic_cast<PointsDumper*>(this))
      ptr->dump_binary(data,size);
    else
      this->dump_binary(data,size);
    break;
  default:
    cRacklet::error("DataDumper do not know how to dump using this OutputFormat");
    break;      
  }
}

/* -------------------------------------------------------------------------- */
inline DataDumper::DataDumper(SpectralModel & mdl) {

  model = & mdl;
}

/* -------------------------------------------------------------------------- */
inline DataDumper::~DataDumper(){

  std::map<std::string,std::ofstream*>::iterator it;
  for (it = timer.begin(); it != timer.end(); ++it){      
    (it->second)->close();
    delete it->second;
  }
  std::map<std::string,Dumper*>::iterator it_d;
  for (it_d = dumper.begin(); it_d != dumper.end(); ++it_d){      
    delete it_d->second;
  }
}
