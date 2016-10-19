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
  E_nor = NULL;
  E_shr = NULL;
  E_fri = NULL;
}

/* -------------------------------------------------------------------------- */
inline DataDumper::~DataDumper(){

  
  if(E_file.is_open())
    E_file.close();
 
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
