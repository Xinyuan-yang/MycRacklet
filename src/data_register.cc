/**
 * @file   data_register.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Aug 17 12:21:46 2016
 *
 * @brief  Implementation of DataRegister class
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
#include "data_register.hh"
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#if defined (_OPENMP)
#include <omp.h>
#endif
/* -------------------------------------------------------------------------- */
// Static members should be defined in a source file
std::string DataRegister::output_dir;
std::string DataRegister::restart_dir;
std::ofstream DataRegister::out_summary;
std::ofstream DataRegister::out_parameters;
std::map<DataFields,DataTypes> DataRegister::datas;
std::map<std::string,Computer*> DataRegister::computers;
std::map<std::string,std::string> DataRegister::variables;
UInt DataRegister::dim;
UInt DataRegister::it;
std::vector<UInt> DataRegister::n_ele;
Real DataRegister::beta;
Real DataRegister::dxmin;
Real DataRegister::ksi;
std::vector<Real> DataRegister::eta;
/* -------------------------------------------------------------------------- */
void DataRegister::data_initialize(const std::string output_folder,
				   const std::string description) {

  output_dir=output_folder;
  restart_dir= "restart_files/";
  extern std::string cR_release_info;

  struct stat info;

  if (stat(output_folder.c_str(),&info) !=0){
    std::stringstream err;
    err << "The output folder does not exists: " << output_folder;
    cRacklet::error(err);
  }
    
  out_summary.open(output_folder+"Simulation_Summary.cra");

  out_summary << " ** This file is automatically generated by cRacklet ** " << std::endl;
  out_summary << std::endl;
  out_summary << " Thank you for using cRacklet !" << std::endl;
  out_summary << " Hereunder is a summary of your last simulation " << std::endl;
  out_summary << std::endl;
  out_summary << "/* -------------------------------------------------------------------------- */ "; 
  out_summary << std::endl;
  out_summary << " * Simulation Description: " << description << std::endl;
  out_summary << std::endl;
#if defined (_OPENMP)
  out_summary <<" * Multi-threaded run with " << omp_get_max_threads() << " threads" << std::endl;
    out_summary << std::endl;
#endif
  out_summary << " * Date and Time:  " << getCurrentDaytime() << std::endl;
  out_summary << " * Sources information: " << std::endl;
  out_summary << cR_release_info << std::endl;
  
  out_parameters.open(output_folder+"Parameters.cra");

#if defined (_OPENMP)
  int nthreads = omp_get_max_threads();
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
#endif

}

/* -------------------------------------------------------------------------- */
void DataRegister::data_finalize() {

  out_summary << std::endl;
  out_summary << "/* -------------------------------------------------------------------------- */ "; 
  out_summary << std::endl;
  out_summary << " Thank you for using cRacklet, hope to see you soon.";
 
  out_summary.close();
  out_parameters.close();

  datas.clear();
  computers.clear();
  variables.clear();

#if defined (_OPENMP)
  fftw_cleanup_threads();
#endif
  fftw_cleanup();
}

/* -------------------------------------------------------------------------- */
std::string DataRegister::getCurrentDaytime() {
  
  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);

}

/* -------------------------------------------------------------------------- */
void DataRegister::registerComputer(std::string computer_name, Computer * computer) {

  std::map<std::string,Computer*>::iterator it = computers.find(computer_name);
  
  if(it == computers.end())
    computers[computer_name] = computer;
  else {
    std::stringstream err;
    err << "Computer named " << computer_name << " already registered !";
    cRacklet::error(err);
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
void DataRegister::registerParameter(std::string name, T value) {

  std::map<std::string,std::string>::iterator it = variables.find(name);

  std::string string_val;
  std::stringstream stream;
  stream << std::fixed;
  stream << std::setprecision(12);
  stream << value;
  stream >> string_val;
  
  if(it == variables.end())
    variables[name] = string_val;
  else {
    std::stringstream err;
    err << "Parameter named " << name << " already registered !";
    cRacklet::error(err);
  }
}

template
void DataRegister::registerParameter(std::string name, Real value);

template
void DataRegister::registerParameter(std::string name, UInt value);

template
void DataRegister::registerParameter(std::string name, std::string value);

template
void DataRegister::registerParameter(std::string name, bool value);

/* -------------------------------------------------------------------------- */
void DataRegister::readInputFile(std::string filename) {
  
  std::ifstream file;
  std::string line;
  file.open(filename);

  if (!file.is_open()){
    std::stringstream err;
    err << "Unable to open input file" << std::endl;
    err << "Check that the file " << filename
	<< " is in the current folder" << std::endl;
    cRacklet::error(err);
  }
    
  while(!file.eof()) { 
    std::getline(file,line);
    std::stringstream sstr(line);
    std::string entry_name;
    std::string parameter;
    sstr >> entry_name;
    if((entry_name == "%")||(line.length()==0))
      continue;
    sstr >> parameter;
    std::cout << entry_name << " has a value of " << parameter << std::endl;
    registerParameter(entry_name, parameter);
  }
}

/* -------------------------------------------------------------------------- */
Computer * DataRegister::getComputer(std::string computer_name) {

  std::map<std::string,Computer*>::iterator it = computers.find(computer_name);

  if(it != computers.end())
    return computers[computer_name];
  else {
    std::stringstream err;
    err << "No registered computer named " << computer_name << " !";
    cRacklet::error(err);
    return NULL;
  }
}

/* -------------------------------------------------------------------------- */
void DataRegister::computeAll(Real time) {

  std::map<std::string,Computer*>::iterator it = computers.begin();

  for (; it!=computers.end(); ++it) {
    (it->second)->compute(time);
  }  
}

/* -------------------------------------------------------------------------- */
void DataRegister::restartComputer(bool pausing,UInt nele_2d) {

  std::ios::openmode mode;

  if(pausing)
    mode=std::ios::out|std::ios::binary;
  else
    mode=std::ios::in|std::ios::binary;
   
  std::fstream file((output_dir+restart_dir+"restart_computers.cra").c_str(),mode);

  if (!file.is_open()&&!pausing)
    cRacklet::error("Unable to open computers restart file");
  
  std::map<std::string,Computer*>::iterator it = computers.begin();

  UInt nb_restarted_computers = 0;
  
  for (; it!=computers.end(); ++it) {
    (it->second)->restart(file,pausing);
    ++nb_restarted_computers;
  }

  if((nele_2d!=0)&&(nb_restarted_computers>0))
    std::cout << "!! WARNING: Computers have been restarted from a 2d simulation !" << std::endl
	      << "!! Be careful with the associated computing domain; "
	      << "2d simulation considers a unit width (i.e. dom_sizez=1[m])"
	      << std::endl; 
  
  file.close();
}

/* -------------------------------------------------------------------------- */
void DataRegister::restart(bool pausing, UInt nele_2d) {

  restartComputer(pausing,nele_2d);
  
  CrackProfile * top_displacements = datas[_top_displacements];
  DataRegister::restartData(top_displacements->getValues(),"restart_top_displacements.cra",pausing, 3*nele_2d);
  CrackProfile * bot_displacements = datas[_bottom_displacements];
  DataRegister::restartData(bot_displacements->getValues(),"restart_bottom_displacements.cra",pausing, 3*nele_2d);

  CrackProfile * top_velocities = datas[_top_velocities];
  DataRegister::restartData(top_velocities->getValues(),"restart_top_velocities.cra",pausing, 3*nele_2d);
  CrackProfile * bot_velocities = datas[_bottom_velocities];
  DataRegister::restartData(bot_velocities->getValues(),"restart_bottom_velocities.cra",pausing, 3*nele_2d);

}
  
/* -------------------------------------------------------------------------- */
template<typename T>
void DataRegister::restartData(std::vector<T> & my_data ,const std::string data_file, bool pausing, UInt nele_2d) {

  std::ios::openmode mode;

  if(pausing)
    mode=std::ios::out|std::ios::binary;
  else
    mode=std::ios::in|std::ios::binary;

  std::string filename = output_dir+restart_dir+data_file;
  
  std::fstream file(filename.c_str(),mode);  

  if (!file.is_open()&&!pausing){
    std::stringstream err;
    err << "Unable to open restart file " << data_file << std::endl;
    cRacklet::error(err);
  }
  
  UInt nb_elements = my_data.size();
  
  if (pausing)
    file.write((char*) &(my_data[0]), nb_elements*sizeof(T));
  else if(nele_2d==0)
    file.read((char*) &(my_data[0]), nb_elements*sizeof(T));
  else
    restartDataFrom2d(my_data,file,nele_2d);
   
  file.close();
}

template
void DataRegister::restartData(std::vector<UInt> & my_data ,const std::string data_file, bool pausing, UInt nele_2d);
template
void DataRegister::restartData(std::vector<Real> & my_data ,const std::string data_file, bool pausing, UInt nele_2d);

/* -------------------------------------------------------------------------- */
template<typename T>
void DataRegister::restartDataFrom2d(std::vector<T> & my_data , std::fstream & file, UInt nele_2d) {

  T * temp_data = new T[nele_2d];
  file.read((char*)(temp_data), nele_2d*sizeof(T));

  for (UInt i = 0; i < my_data.size(); i+=nele_2d) {
    memcpy(&(my_data[i]),temp_data,nele_2d*sizeof(T));
  }
  delete[] temp_data;
}

template
void DataRegister::restartDataFrom2d(std::vector<UInt> & my_data , std::fstream & file, UInt nele_2d);
template
void DataRegister::restartDataFrom2d(std::vector<Real> & my_data , std::fstream & file, UInt nele_2d);

/* -------------------------------------------------------------------------- */
UInt DataRegister::getCrackTipPosition(UInt x_start, UInt x_end, UInt z_pos) {

  UInt nb_x = n_ele[0];
  
  UInt x_tip=x_start;

  const std::vector<UInt> * ind_crack = readData(_id_crack);

  while (((*ind_crack)[x_tip+z_pos*nb_x]==2)&&(x_tip<x_end))
    ++x_tip;
  
  return x_tip;
}

/* -------------------------------------------------------------------------- */
UInt DataRegister::getCohesiveTipPosition(UInt x_start, UInt x_end, UInt z_pos) {

  UInt nb_x = n_ele[0];
  
  UInt x_tip=x_start;

  const std::vector<UInt> * ind_crack = readData(_id_crack);

  while (((*ind_crack)[x_tip+z_pos*nb_x]==1)&&(x_tip<x_end))
    ++x_tip;
  
  return x_tip;
}

/* -------------------------------------------------------------------------- */
void Integrator::integrate(Real time){

  Real dt = time - t_old;
  
  I += 0.5*dt*(I_dot_old + I_dot);
  I_dot_old = I_dot;
  I_dot = 0;
  t_old = time;
}

/* -------------------------------------------------------------------------- */
void SurfingIntegrator::updateIntegrationPoints() {
  
  UInt crack_tip = DataRegister::getCrackTipPosition(crack_start,crack_end);

  UInt start = std::max(0,(int)(crack_tip-integ_width));
  UInt end = std::max(0,(int)(crack_tip+integ_width));

  std::vector<UInt> new_integ_points(end-start);

  for (UInt i = 0; i < (end-start); ++i) {
    new_integ_points[i] = start+i;
  }
  this->index = new_integ_points;  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getTopDisplacements(){
  
  const CrackProfile * u_top = DataRegister::readData(_top_displacements);

  return u_top;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getBotDisplacements(){

  const CrackProfile * u_bot = DataRegister::readData(_bottom_displacements);
  return u_bot;

}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getShearDisplacementJumps(){
  
  const CrackProfile * shear_u_jump = DataRegister::readData(_shear_displacement_jumps);
  return shear_u_jump;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getNormalDisplacementJumps(){
  
  const CrackProfile * normal_u_jump = DataRegister::readData(_normal_displacement_jumps);

  return normal_u_jump;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getTopVelocities(){
  
  const CrackProfile * v_top = DataRegister::readData(_top_velocities);

  return v_top;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getBotVelocities(){

  const CrackProfile * v_bot = DataRegister::readData(_bottom_velocities);
  return v_bot;

}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getShearVelocityJumps(){
  
  const CrackProfile * shear_v_jump = DataRegister::readData(_shear_velocity_jumps);
  return shear_v_jump;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getNormalVelocityJumps(){
  
  const CrackProfile * normal_v_jump = DataRegister::readData(_normal_velocity_jumps);

  return normal_v_jump;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getInterfaceTractions(){
  
  const CrackProfile * itf_tractions = DataRegister::readData(_interface_tractions);

  return itf_tractions;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getTopDynamicStresses(){
  
  const CrackProfile * top_dyn_stresses = DataRegister::readData(_top_dynamic_stress);

  return top_dyn_stresses;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getBotDynamicStresses(){
  
  const CrackProfile * bot_dyn_stresses = DataRegister::readData(_bottom_dynamic_stress);

  return bot_dyn_stresses;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getTopLoading(){
  
  const CrackProfile * top_loading = DataRegister::readData(_top_loading);

  return top_loading;
  
}

/* -------------------------------------------------------------------------- */
const CrackProfile * DataRegister::getBotLoading(){
  const CrackProfile * bot_loading = DataRegister::readData(_bottom_loading);
  return bot_loading;
}
