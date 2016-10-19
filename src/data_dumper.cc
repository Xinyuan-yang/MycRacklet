/* -------------------------------------------------------------------------- */
#include "data_dumper.hh"
#include <time.h>
/* -------------------------------------------------------------------------- */
void Dumper::dump() {
  
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
  case _top_loading:
  case _bottom_loading:
    dump<CrackProfile>();
    break;

  case _normal_strength:
  case _shear_strength:
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
void PointsDumper::dump() {

  for (UInt i = 0; i < points.size(); ++i) {
    this->start = points[i];
    file << points[i] <<  " ";    
    for (UInt f = 0; f < fields.size(); ++f) {
      this->field_name = fields[f];
      Dumper::dump();
    }
  }
  file << std::endl;
}

/* -------------------------------------------------------------------------- */
void DataDumper::dump(std::string filename) {
  dumper[filename]->dump();
  dumper[filename]->endl();
  printTime(filename);
  
}

/* -------------------------------------------------------------------------- */
void DataDumper::dumpAll() {

  std::map<std::string,Dumper*>::iterator it;
  for (it = dumper.begin(); it != dumper.end(); ++it){      
    (it->second)->dump();
    (it->second)->endl();
    printTime(it->first);
  }
}

/* -------------------------------------------------------------------------- */
void DataDumper::initEnergetics(const std::string filename) {

  E_nor = & model->getNormalDissipatedEnergy();
  E_shr = & model->getShearDissipatedEnergy();
  E_fri = & model->getFrictionalEnergy();
  Eq = & model->getRadiatedEnergy();

  E_file.open((DataRegister::output_dir+filename).c_str());
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "ENERGY OUTPUT FILE: " << std::endl << filename
 			    << std::endl;
  
}

/* -------------------------------------------------------------------------- */
void DataDumper::initDumper(const std::string filename, DataFields type, Real nb_data_ratio, UInt stride, UInt start,OutputFormat format) {

  dumper[filename] = new Dumper(type, nb_data_ratio,stride,start);
  dumper[filename]->init(DataRegister::output_dir+filename,format);
  initTimer(filename);
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "FIELD OUTPUT FILE: " << std::endl << filename
			    << " (dumping ratio: " << nb_data_ratio
			    << ", stride: " << stride
			    << ", start: " << start << ")"
			    << std::endl;

}

/* -------------------------------------------------------------------------- */
void DataDumper::initVectorDumper(const std::string filename, DataFields type, UInt dimension_to_dump,
			Real ratio_of_nele, UInt stride, UInt start, OutputFormat format) {

  this->initDumper(filename,type,ratio_of_nele/3.0,stride*3,start+dimension_to_dump,format);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsDumper(const std::string filename, std::vector<DataFields> & fields,
				  std::vector<UInt> & points_to_dump) {
  UInt total_n_ele = model->getNbElements()[0]*model->getNbElements()[1];
  dumper[filename] = new PointsDumper(fields, points_to_dump, total_n_ele);
  dumper[filename]->init(DataRegister::output_dir+filename);
  initTimer(filename);
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "POINTWISE OUTPUT FILE: " << std::endl << filename << std::endl;

}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsDumper(const std::string filename, std::vector<UInt> & points_to_dump) {
  this->initPointsDumper(filename, standard_fields_history, points_to_dump);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsDumper(const std::string filename, UInt start, UInt end, UInt nb_obs_points) {
  UInt point_spacing = (end-start)/nb_obs_points;
  UInt pos = start+(UInt)(0.5*point_spacing);
  
  std::vector<UInt> point_his(nb_obs_points);

  for (UInt p = 0; p < nb_obs_points; ++p) {
    point_his[p] = pos;
    pos += point_spacing; 
  }  
}
/* -------------------------------------------------------------------------- */
void DataDumper::initTimer(std::string filename) {

  std::string timer_filename = "Timer_"+filename;
  
  timer[filename]=new std::ofstream;
  timer[filename]->open(DataRegister::output_dir+timer_filename);
}

/* -------------------------------------------------------------------------- */
void DataDumper::printEnergetics() {
  
  E_file << std::scientific << std::setprecision(9)
	 << model->getTime() <<" "<< (*E_shr)[0].E+(*E_nor)[0].E <<" " 
	 <<(*E_shr)[0].E_dot_old+(*E_nor)[0].E_dot_old <<" "<< (*E_fri)[0].E <<" "
	 << (*E_fri)[0].E_dot_old <<" "<< (*E_shr)[1].E+(*E_nor)[1].E <<" "
	 << (*E_shr)[1].E_dot_old+(*E_nor)[1].E_dot_old <<" "<< (*E_fri)[1].E <<" "
	 << (*E_fri)[1].E_dot_old <<" "<< (*E_nor)[0].E <<" "<< (*E_nor)[0].E_dot_old <<" "
	 << (*E_shr)[0].E <<" "<< (*E_shr)[0].E_dot_old <<" "<<	(*E_nor)[1].E <<" "
	 << (*E_nor)[1].E_dot_old <<" "<< (*E_shr)[1].E <<" "<< (*E_shr)[1].E_dot_old <<" "
	 << Eq->E << " " << Eq->E_dot_old << " " << std::endl;
}

/* -------------------------------------------------------------------------- */
void DataDumper::printTime(std::string filename) {
 
  *(timer[filename]) << std::scientific << std::setprecision(9)
		     << model->getTime() << std::endl;
}

