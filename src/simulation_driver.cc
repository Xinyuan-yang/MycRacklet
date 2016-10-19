/* -------------------------------------------------------------------------- */
#include "simulation_driver.hh"
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
SimulationDriver::SimulationDriver(SpectralModel & model_to_drive,
				   Real target_speed, Real crack_start)
  : model(model_to_drive) { 

  this->target_crack_speed=target_speed;
  double new_beta = this->adjustStableTimeStep();
  nb_t_by_elem=1/(new_beta*target_crack_speed);
  model.initModel(new_beta);
  std::vector<Real> dx = model.getElementSize();
  this->x_crack_start = (UInt)(crack_start/dx[0]);  
}

/* -------------------------------------------------------------------------- */
void SimulationDriver::initConstantLoading(Real initial_loading, Real psi, Real phi) {
  this->target_crack_speed=0.;
  model.setLoadingCase(initial_loading, psi, phi);
  model.updateLoads();
  model.computeInitialVelocities();
}
/* -------------------------------------------------------------------------- */
void SimulationDriver::initConstantSpeed(Real initial_loading, Real psi, Real phi,
					 Real average_crit_stress,
					 LoadControlType load_control) {

  this->lc_type=load_control;
  std::vector<UInt> * id_crack = datas[_id_crack];
  model.setLoadingCase(initial_loading, psi, phi);
  const std::vector<Real> & uniform_loading = model.getUniformLoading();
  new_loading.resize(model.getDim());
  for (UInt i = 0; i < model.getDim(); ++i) {
    new_loading[i] = uniform_loading[i];
  }
  
  this->av_crit_stress=average_crit_stress;

  std::vector<UInt> n_ele = model.getNbElements();
  std::vector<Real> dx = model.getElementSize();
  UInt x_tip = model.getCrackTipPosition(this->x_crack_start);

  std::cout << "Crack position initially detected at position: " << x_tip*dx[0] << std::endl;
  //constant crack propagation speed loading init.
  at.resize(2);
  at[0] = (x_tip+5)*dx[0];
  this->c_s_initiation = true;
  countcfs=1;
  atc.resize(2);
  atc[0]=0;
  
  switch (lc_type) {
  case _time_control:
    setLoadingCase<_time_control>();
    break;
  case _space_control:
    setLoadingCase<_space_control>();
    break;
  default:
    cRacklet::error("Please specify a valid crack speed control algorithm.");
  }
  model.computeInitialVelocities();
}

/* -------------------------------------------------------------------------- */
Real SimulationDriver::initLoadingFromFile(std::string loading_file,
					   LoadControlType load_control) {

  Real nb_steps_file;
  this->lc_type=load_control;
  std::ifstream file (loading_file,std::ios::in);
  if (!file.is_open()){
    std::stringstream err;
    err << "Unable to open file" << std::endl;
    err << "Check that the file" << loading_file 
	<<" is in the current folder" << std::endl;
    cRacklet::error(err);
  }
  
  switch (lc_type) {
  case _time_control:
    nb_steps_file = readSetLoadingCase<_time_control>(file);
    break;
  case _space_control:
    nb_steps_file = readSetLoadingCase<_space_control>(file);
    break;
  default:
    cRacklet::error("Please specify a valid .");
  }
  
  runReadingStep();
  this->target_crack_speed=-1;

  model.computeInitialVelocities();
    
  return nb_steps_file;
}
/* -------------------------------------------------------------------------- */
void SimulationDriver::checkForTargetSpeed(std::ifstream & file) {

  Real v_obj;
  file >> v_obj;

  if (v_obj==target_crack_speed)  
    std::cout << "Loading tailored to fix crack speed at " << v_obj << "c_s has just been read." << std::endl;
  else {
    std::stringstream error;
    error << "!!! A problem occured with your loading file! Read target speed " << v_obj 
	  << " instead of " << target_crack_speed;
    cRacklet::error(error);
  }
}
/* -------------------------------------------------------------------------- */
Real SimulationDriver::adjustStableTimeStep() {

  UInt nb_it=std::ceil(1/(target_crack_speed*CS_DT_OVER_DX));
  Real new_beta = 1.0/(nb_it*target_crack_speed);
  return new_beta;
}

void SimulationDriver::priorUpdates() {

  model.updateDisplacements();
  model.updateMaterialProp();
  
}
/* -------------------------------------------------------------------------- */
void SimulationDriver::solveTimeStep() {

  model.fftOnDisplacements();
  model.computeStress();
  model.computeVelocities();
  model.computeEnergy();
  model.increaseTimeStep();
}

/* -------------------------------------------------------------------------- */
UInt SimulationDriver::solveStep() {

  UInt ret=0;

  priorUpdates();
  
  if(this->target_crack_speed) {
    if(this->target_crack_speed==-1.)
      ret = runReadingStep();
    else
      ret = runWritingStep();
  }
  
  solveTimeStep();

  return ret;  
}

/* -------------------------------------------------------------------------- */
UInt SimulationDriver::runReadingStep() {

  UInt x_tip = model.getCrackTipPosition(x_crack_start);
  UInt dim = model.getDim();
  UInt it = model.getCurrentTimeStep();
  
  switch (lc_type) {
  case _time_control:
    if (it > ctrl_loading.size()/3)
      model.updateLoads(&(ctrl_loading[dim*(ctrl_loading.size()/3-1)]));
    else
      model.updateLoads(&(ctrl_loading[dim*(it-1)]));
    break;
  case _space_control:
    model.updateLoads(&(ctrl_loading[dim*(x_tip-x_crack_start)]));
    break;
  default:
    cRacklet::error("Ill-initialized crack speed control algorithm. No type valid specified.");
  }
 
  return x_tip;
}

/* -------------------------------------------------------------------------- */
UInt SimulationDriver::runWritingStep() {

  UInt x_tip;
  bool has_progressed = cstCrackFrontSpeed(x_tip);
  UInt it = model.getCurrentTimeStep();
  
  switch (lc_type) {
  case _time_control:
    if(model.getNbTimeSteps()>0) {
      for (UInt j = 0; j < model.getDim(); ++j) {
	ctrl_loading[3*(it-1)+j]=new_loading[j];
      }
    }
    else {
      for (UInt j = 0; j < model.getDim(); ++j) {
	ctrl_loading.push_back(new_loading[j]);
      }
    }
    break;
  case _space_control:
    if (has_progressed) {
      for (UInt j=0; j < model.getDim(); ++j) {
	ctrl_loading[3*x_tip+j] = new_loading[j];
      }
    }
    break;
  default:
    cRacklet::error("Ill-initialized crack speed control algorithm. No type valid specified.");
  }
  model.updateLoads(&(new_loading[0]));
  return x_tip;
}

/* -------------------------------------------------------------------------- */
bool SimulationDriver::cstCrackFrontSpeed(UInt & x_tip){

  bool has_progressed=false;

  std::vector<Real> dx = model.getElementSize();
  x_tip = model.getCrackTipPosition(this->x_crack_start);
  at[1] = (x_tip-0.5)*dx[0];

  std::vector<Real> strr;
  strr.resize(3);
  strr[0]=new_loading[0]/av_crit_stress;
  strr[1]=new_loading[1]/av_crit_stress;
  strr[2]=new_loading[2]/av_crit_stress;
    
  //check if the cohesive zone has propagated 
  if (at[0]<at[1]) {
      //check if crack has reached target speed during initiation
      if((nb_t_by_elem>=countcfs)||(!c_s_initiation)) {
	if(c_s_initiation) {
	  std::cout << "Initiation phase just ends -> time: " << model.getTime() << std::endl;
	  c_s_initiation = false;
	}
	atc[1]=countcfs/nb_t_by_elem;
	Real rapp;
	if (atc[0]==0||(atc[1]<1.01&&atc[1]>0.99)) {
	  rapp=atc[1];
	}
	else {rapp=(atc[1]+atc[0])/2;atc[1]=rapp;}
      
	for (UInt j=0; j < model.getDim(); ++j) {
	  if (strr[j]*rapp>=0.525) {
	    new_loading[j]=0.525*av_crit_stress;
	  }
	  else{
	    new_loading[j]=new_loading[j]*rapp;
	  }
	}
	has_progressed=true;
      }
      countcfs=1;
      at[0]=at[1];
    }
    else{countcfs++;}
  
  atc[0]=atc[1];
  
  return has_progressed;
}

/* -------------------------------------------------------------------------- */
void SimulationDriver::writeLoading(std::string load_file ) {

  std::vector<UInt> nb_ele = model.getNbElements();
  UInt nb_something;
  switch (lc_type) {
  case _time_control:
    nb_something = model.getNbTimeSteps();
    break;
  case _space_control:
    nb_something = nb_ele[0];
    break;
  }
  std::ofstream os;
  os.open(load_file);
  os.precision(16);
  os << std::scientific;
  os << nb_something << std::endl;
  os << this->target_crack_speed << std::endl;
  for (UInt i=0; i < ctrl_loading.size(); ++i) {
    os <<ctrl_loading[i]<<std::endl;
  }
  os.close();
}

/* -------------------------------------------------------------------------- */
void SimulationDriver::launchCrack(Real crack_start, Real launched_size,
				   Real v_init) {

  std::vector<Real> dx = model.getElementSize();
  
  UInt l_end = (UInt)(launched_size/dx[0]);
  UInt x_start = (UInt)(crack_start/dx[0]);
  UInt every_t = (UInt)(1/(model.getBeta()*v_init));
  UInt x_tip = model.getCrackTipPosition(x_start);

  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  
  std::cout << "Rupture is currently artificially triggered with speed " << v_init
	    << "*c_s" << std::endl; 
  
  while ((x_tip-x_start)<l_end) {

    priorUpdates();
    
    x_tip = model.getCrackTipPosition(x_start);

    if (model.getCurrentTimeStep()%every_t==0) {
      (*nor_strength)[x_tip]=0;
      (*shr_strength)[x_tip]=0;
    }
    solveTimeStep();
  }
  std::cout << "Initiation just ends..." << std::endl; 
}
