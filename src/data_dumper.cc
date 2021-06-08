/**
 * @file   data_dumper.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Sep  3 17:04:05 2014
 *
 * @brief  Implementation of the DataDumper class and Dumper class
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
#include "data_dumper.hh"
#include <time.h>
/* -------------------------------------------------------------------------- */
void PointsDumper::dump() {

  for (UInt i = 0; i < points.size(); ++i) {
    this->start = points[i];
    for (UInt f = 0; f < fields.size(); ++f) {
      this->field_name = fields[f];
      Dumper::dump();
    }
  }
}

/* -------------------------------------------------------------------------- */
void IntegratorsDumper::dump() {

  std::vector<Integrator*>::iterator it;
  for (it = integrators.begin(); it != integrators.end(); ++it){      
    file << std::scientific << std::setprecision(9)
	 << (*it)->getIntegration()
	 << " "
	 << (*it)->getRate()
	 << " ";
  }
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
void DataDumper::initDumper(const std::string filename, DataFields type, Real nb_data_ratio, UInt stride, UInt start,OutputFormat format) {

  dumper[filename] = new Dumper(type, nb_data_ratio,stride,start);
  dumper[filename]->init(DataRegister::output_dir+filename,format);
  initTimer(filename);
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "FIELD OUTPUT FILE OF " << datafields_name[type] << std::endl
			    << filename
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
void DataDumper::initPointsDumper(const std::string filename, std::vector<DataFields> fields,
				  std::vector<UInt> points_to_dump, OutputFormat format) {
  UInt total_n_ele = model->getNbElements()[0]*model->getNbElements()[1];
  dumper[filename] = new PointsDumper(fields, points_to_dump, total_n_ele);
  dumper[filename]->init(DataRegister::output_dir+filename, format);
  initTimer(filename);
  
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "POINTWISE OUTPUT FILE OF: ";

  std::vector<DataFields>::iterator it;
  for(it=fields.begin(); it!=fields.end(); ++it)
    DataRegister::out_summary << datafields_name[*it] << "; ";

  DataRegister::out_summary << std::endl << filename << std::endl;

  //Printing details of observation points positions
  std::string info_file = DataRegister::output_dir+"Info_"+filename;
  std::ofstream file(info_file);
  std::vector<UInt>::iterator itp;
  UInt nele_x = model->getNbElements()[0];
  file << points_to_dump.size() << std::endl;
  for(itp=points_to_dump.begin(); itp!=points_to_dump.end(); ++itp) {
    UInt iz = *itp/nele_x;
    UInt ix = *itp%nele_x;
    file << ix << " " << iz << std::endl;
  }
  file.close();
}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsDumper(const std::string filename, std::vector<UInt> points_to_dump,
				  OutputFormat format) {
  this->initPointsDumper(filename, standard_fields_history, points_to_dump);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsDumper(const std::string filename, UInt start, UInt end, UInt nb_obs_points,
				  OutputFormat format) {
  UInt point_spacing = (end-start)/nb_obs_points;
  UInt pos = start+(UInt)(0.5*point_spacing);
  
  std::vector<UInt> point_his(nb_obs_points);

  for (UInt p = 0; p < nb_obs_points; ++p) {
    point_his[p] = pos;
    pos += point_spacing; 
  }
  this->initPointsDumper(filename, point_his);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initIntegratorsDumper(const std::string filename,
				       std::vector<UInt> integration_domain,
				       std::vector<IntegratorTypes> inte_types,
				       std::vector<std::string> integrator_names,
				       OutputFormat format) {
 
  Real time = model->getTime();
  std::vector<Real> dx = model->getElementSize();
  IntegratorsDumper * ptr = new IntegratorsDumper(inte_types, integrator_names, integration_domain,
						  time, dx[0]*dx[1]);
  dumper[filename] = ptr;
  dumper[filename]->init(DataRegister::output_dir+filename, format);
  initTimer(filename);
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "INTEGRATORS OUTPUT FILE: " << std::endl << filename << std::endl;
  for (UInt i = 0; i < integrator_names.size(); ++i) {
    DataRegister::out_summary << "* " << i << ": " << integrator_names[i] << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
void DataDumper::initIntegratorsDumper(const std::string filename,
				       std::vector<Real> start_corner, std::vector<Real> end_corner,
				       OutputFormat format) {
  
  std::vector<UInt> n_ele = model->getNbElements();
  std::vector<Real> dx = model->getElementSize();
  
  std::vector<UInt> start(2);
  std::vector<UInt> end(2);

  for (UInt i = 0; i < 2; ++i) {
    start[i] = (UInt)(start_corner[i]/dx[i]);
    end[i] = (UInt)(end_corner[i]/dx[i]);
  }

  if(n_ele[1]==1){
    start[1]=0;
    end[1]=1;
  }
  
  std::vector<UInt> int_points;
  
  for (UInt ix = start[0]; ix < std::min(n_ele[0],end[0]); ++ix) {
    for (UInt iz = start[1]; iz < std::min(n_ele[1],end[1]); ++iz) {
      int_points.push_back(ix+iz*n_ele[0]);
    }
  }

  this->initIntegratorsDumper(filename,int_points,
			      standard_energy_integrators,standard_energy_names,format);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initIntegratorsDumper(const std::string filename, OutputFormat format) {

  std::vector<UInt> n_ele = model->getNbElements();
  std::vector<Real> dx = model->getElementSize();

  this->initIntegratorsDumper(filename,{0.,0.},{n_ele[0]*dx[0],n_ele[1]*dx[1]},format);  
}

/* -------------------------------------------------------------------------- */
void DataDumper::initSurfingIntegratorsDumper(const std::string filename,  UInt integration_width,
					      UInt crack_start, UInt crack_end,
					      std::vector<IntegratorTypes> inte_types,
					      std::vector<std::string> integrator_names, OutputFormat format){
  
  Real time = model->getTime();
  std::vector<Real> dx = model->getElementSize();
  IntegratorsDumper * ptr = new IntegratorsDumper(inte_types, integrator_names, integration_width,
						  time, dx[0]*dx[1],crack_start,crack_end);
  dumper[filename] = ptr;
  dumper[filename]->init(DataRegister::output_dir+filename, format);
  initTimer(filename);
  DataRegister::out_summary << "/* -------------------------------------------------------------------------- */"; 
  DataRegister::out_summary << std::endl;
  DataRegister::out_summary << "SURFING INTEGRATORS OUTPUT FILE: " << std::endl << filename << std::endl;
  DataRegister::out_summary << "Integration width: " << std::endl << integration_width << std::endl;
  for (UInt i = 0; i < integrator_names.size(); ++i) {
    DataRegister::out_summary << "* " << i << ": " << integrator_names[i] << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
void DataDumper::initTimer(std::string filename) {

  std::string timer_filename = "Timer_"+filename;
  
  timer[filename]=new std::ofstream;
  timer[filename]->open(DataRegister::output_dir+timer_filename);
}

/* -------------------------------------------------------------------------- */
void DataDumper::printTime(std::string filename) {
 
  *(timer[filename]) << std::scientific << std::setprecision(9)
		     << model->getTime() << std::endl;
}

