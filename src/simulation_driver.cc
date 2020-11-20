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
  model.initInterfaceFields();
}
/* -------------------------------------------------------------------------- */
void SimulationDriver::initConstantSpeed(Real initial_loading, Real psi, Real phi,
					 Real average_max_stress, Real spont_crack_length,
					 LoadControlType load_control,
					 Real load_upper_bound, Real griffith_length) {

  this->lc_type=load_control;
  this->spont_crack=spont_crack_length;
  this->max_load=load_upper_bound;
  this->min_load=griffith_length*initial_loading*initial_loading;
  model.setLoadingCase(initial_loading, psi, phi);
  const std::vector<Real> & uniform_loading = model.getUniformLoading();
  new_loading.resize(model.getDim());
  for (UInt i = 0; i < model.getDim(); ++i) {
    new_loading[i] = uniform_loading[i];
  }
  
  this->av_max_stress=average_max_stress;

  std::vector<UInt> n_ele = model.getNbElements();
  std::vector<Real> dx = model.getElementSize();

  //constant crack propagation speed loading init.
  at.resize(2);
  at[0] = 0;
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
  model.updateLoads();
  model.initInterfaceFields();
}

/* -------------------------------------------------------------------------- */
Real SimulationDriver::initLoadingFromFile(std::string loading_file,
					   LoadControlType load_control,
					   Real initial_loading, Real psi, Real phi) {

  Real nb_steps_file;
  this->lc_type=load_control;
  std::ifstream file (loading_file,std::ios::in);
  if (!file.is_open()){
    std::stringstream err;
    err << "Unable to open file" << std::endl;
    err << "Check that the file " << loading_file 
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
 
  this->target_crack_speed=-1;

  if(initial_loading != 0.) {
    model.setLoadingCase(initial_loading, psi, phi);
    model.updateLoads();
  }
  else {
    runReadingStep();
  }

  model.initInterfaceFields();  
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

/* -------------------------------------------------------------------------- */
void SimulationDriver::solveTimeStep() {

  model.updateDisplacements(); 
  model.fftOnDisplacements();
  model.computeStress();
  model.computeInterfaceFields();
  model.increaseTimeStep();

}

/* -------------------------------------------------------------------------- */
UInt SimulationDriver::solveStep() {

  UInt ret=0;
   
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

  UInt x_tip = model.getCrackTipPosition(x_crack_start,model.getNbElements()[0]);
  UInt dim = model.getDim();
  UInt it = model.getCurrentTimeStep();
  
  switch (lc_type) {
  case _time_control:
    if(reading_time==-1)
      reading_time=it-1;
    if ((it-reading_time) > ctrl_loading.size()/3)
      model.updateLoads(&(ctrl_loading[dim*(ctrl_loading.size()/3-1)]));
    else
      model.updateLoads(&(ctrl_loading[dim*(it-reading_time-1)]));
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
  
  switch (lc_type) {
  case _time_control:
    for (UInt j = 0; j < model.getDim(); ++j) {
      ctrl_loading.push_back(new_loading[j]);
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
  x_tip = model.getCrackTipPosition(this->x_crack_start,model.getNbElements()[0]);
  at[1] = (x_tip-0.5)*dx[0];
  
  if(at[0]==0) //First algorithm step
    at[0]=at[1];
   
  //check if the cohesive zone has propagated 
  if (at[0]<at[1]) {
    //check if crack has reached target speed during initiation
    if(((nb_t_by_elem>=countcfs)&&(at[1]>spont_crack))||(!c_s_initiation)) {
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

      Real res_tau=0.;
      for (UInt j=0; j < model.getDim(); ++j) {
	res_tau += new_loading[j]*new_loading[j]; 
      }
      res_tau = sqrt(res_tau);
      rapp = std::min(max_load*av_max_stress/res_tau,rapp);
      rapp = std::max(sqrt(min_load/(x_tip*dx[0]))/res_tau,rapp);
      
      for (UInt j=0; j < model.getDim(); ++j) {
	new_loading[j]=new_loading[j]*rapp;
      }

      has_progressed=true;
    }
    countcfs=1;
    at[0]=at[1];
    atc[0]=atc[1];
  }
  else{countcfs++;}

  return has_progressed;
}

/* -------------------------------------------------------------------------- */
void SimulationDriver::writeLoading(std::string load_file ) {

  std::vector<UInt> nb_ele = model.getNbElements();
  UInt nb_something;
  switch (lc_type) {
  case _time_control:
    nb_something = ctrl_loading.size()/3;
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
				   Real v_init, bool one_side_propagation) {

  std::vector<Real> dx = model.getElementSize();
  std::vector<UInt> nb_elements = model.getNbElements();
  UInt l_end = (UInt)(launched_size/dx[0]);
  UInt x_start = (UInt)(crack_start/dx[0]);
  UInt every_t = (UInt)(1/(model.getBeta()*v_init));
  UInt x_tip = model.getCrackTipPosition(x_start,model.getNbElements()[0]);
  UInt x_tip_prev=0;
  
  std::vector<Real> * nor_strength = datas[_normal_strength];
  std::vector<Real> * shr_strength = datas[_shear_strength];
  std::vector<UInt> * ind_crack = datas[_id_crack];
  
  std::cout << "Rupture is currently artificially triggered with speed " << v_init
	    << "*c_s" << std::endl; 

  UInt growth_factor;

  if(one_side_propagation)
    growth_factor=1;
  else
    growth_factor=2;
  
  l_end /= growth_factor;

  while ((x_tip-x_start)<l_end) {
    
    x_tip = model.getCrackTipPosition(x_start,model.getNbElements()[0]);
    
    if (model.getCurrentTimeStep()%every_t==0) {

      if (((x_tip-x_tip_prev)>2)&&(x_tip_prev!=0.)) {
	break;
      }

      std::cout << " Crack position is now at " << growth_factor*(x_tip-x_start)*dx[0] << std::endl;
      x_tip_prev = model.getCrackTipPosition(x_start,model.getNbElements()[0]);
      
      for (UInt z = 0; z < nb_elements[1]; ++z) {
	(*nor_strength)[x_tip+z*nb_elements[0]]=0.;
	(*shr_strength)[x_tip+z*nb_elements[0]]=0.;
 	(*ind_crack)[x_tip+z*nb_elements[0]]=2;

	if(!one_side_propagation) {
	  (*nor_strength)[2*x_start-x_tip+z*nb_elements[0]]=0.;
	  (*shr_strength)[2*x_start-x_tip+z*nb_elements[0]]=0.;
	  (*ind_crack)[2*x_start-x_tip+z*nb_elements[0]]=2;
	}
      }
    }
    solveTimeStep();
  }
  std::cout << "Initiation just ends..." << std::endl; 
}
