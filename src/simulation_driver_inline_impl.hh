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
