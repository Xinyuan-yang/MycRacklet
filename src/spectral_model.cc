/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <gsl/gsl_sf_bessel.h>
#include <complex>
#include <fftw3.h>
#include <math.h>
#include <sys/stat.h>
#include <algorithm>
/* -------------------------------------------------------------------------- */
InterfaceFields::InterfaceFields(const std::vector<CrackProfile> * in_fields, UInt dimension) {

  this->init(&((*in_fields)[0]), &((*in_fields)[1]), dimension);
}

/* -------------------------------------------------------------------------- */
InterfaceFields::InterfaceFields(const CrackProfile* in_fields_top,
				 const CrackProfile* in_fields_bot, UInt dimension) {
  
  
  this->init(in_fields_top, in_fields_bot, dimension);
}

/* -------------------------------------------------------------------------- */
void InterfaceFields::init(const CrackProfile* in_fields_top,
			   const CrackProfile* in_fields_bot, UInt dimension) {

  fields.resize(2);
  fields[0] = in_fields_top;
  fields[1] = in_fields_bot;  
  dim = dimension;
  delta_fields.resize(2);
}

/* -------------------------------------------------------------------------- */
void InterfaceFields::computeJumpFields() {

  fields_jump = *(fields[0]) - *(fields[1]);

  CrackProfile delta_field_1 = fields_jump.getStridedPart(dim,0);
  CrackProfile delta_field_3 = fields_jump.getStridedPart(dim,2);

  delta_fields[0] = delta_field_1*delta_field_1 + delta_field_3*delta_field_3;
  delta_fields[0].squareRoot();

  delta_fields[1] = fields_jump.getStridedPart(dim,1);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::initModel(Real reset_beta, bool blank) {

  dim = 3;
  
  //elements
  total_n_ele = n_ele[0]*n_ele[1];

  dx[0] = X[0]/(Real)(n_ele[0]);
  q0[0] = 2*M_PI/(dx[0]*(Real)n_ele[0]);
  nele_fft[0]=n_ele[0]/2+1;
  nele_fft[1]=n_ele[1];
  
  //for 2D 
  if (n_ele[1]==1) {dxmin=dx[0]; dx[1] = 1.; q0[1]=0.;}
  else if (n_ele[0]==1)
    cRacklet::error("1D Interface are build along x-axis by convention");
  //3D
  else{
    dx[1] = X[1]/(Real)(n_ele[1]);
    q0[1] = 2*M_PI/(dx[1]*(Real)n_ele[1]);
    dxmin=std::min(dx[0],dx[1]);
  }

  if(reset_beta>0.)
    resetStableTimeStep(reset_beta);
  else
    beta=CS_DT_OVER_DX;

  //modes
  total_nele_fft=nele_fft[0]*nele_fft[1];

  it = 0;

  nb_kernels = 4;

  displacements.resize(2);
  velocities.resize(2);
  stresses.resize(2);
  loads.resize(2);
  eta.resize(2);

  loading_ratio.resize(total_n_ele,1.0);

  if(!blank) {
    for (UInt i = 0; i < 2; ++i) {
      displacements[i].SetGridSize(n_ele, dim);
      velocities[i].SetGridSize(n_ele, dim);
      stresses[i].SetGridSize(n_ele, dim);
      loads[i].SetGridSize(n_ele, dim);
      eta[i] = sqrt(2*(1-nu[i])/(1-2*nu[i]));

      displacements[i].initFFT(true,dim);
      stresses[i].initFFT(false,dim);
    }

    intfc_trac.SetGridSize(n_ele, dim);
  }
  
  displ_jump = new InterfaceFields(&displacements, dim);
  veloc_jump = new InterfaceFields(&velocities, dim);
  
  /* -------------------------------------------------------------------------- */
  //set the modes q0, k and m
  if(interface_dim==2)
    initFrequency<2>();
  else if(interface_dim==1)
    initFrequency<1>();

  if(!blank) {
    h11u1 = new Real[2*(total_nele_fft-1)];
    h22u2 = new Real[2*(total_nele_fft-1)];
    h33u3 = new Real[2*(total_nele_fft-1)];
    h12u1 = new Real[2*(total_nele_fft-1)];
    h12u2 = new Real[2*(total_nele_fft-1)];

    U_top = new Real[2*(total_nele_fft-1)*dim];
    U_bot = new Real[2*(total_nele_fft-1)*dim];
    F_k = new Real[2*total_nele_fft*dim];
  }
  else
    h11u1 = NULL;
  
  this->nb_kernels = 4;
  
  initConvolutionManagers(blank);

  if(!blank) {
    registerModelFields();
    printSelf();
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::registerModelFields() {

  this->registerData(_top_displacements, &(displacements[0]));
  this->registerData(_bottom_displacements, &(displacements[1]));
  this->registerData(_shear_displacement_jumps, &(displ_jump->delta_fields[0]));
  this->registerData(_normal_displacement_jumps, &(displ_jump->delta_fields[1])); 
  this->registerData(_top_velocities, &(velocities[0]));
  this->registerData(_bottom_velocities, &(velocities[1]));
  this->registerData(_shear_velocity_jumps, &(veloc_jump->delta_fields[0]));
  this->registerData(_normal_velocity_jumps, &(veloc_jump->delta_fields[1]));
  this->registerData(_top_dynamic_stress, &(stresses[0]));
  this->registerData(_bottom_dynamic_stress, &(stresses[1]));
  this->registerData(_interface_tractions, &intfc_trac);
  this->registerData(_top_loading, &(loads[0]));
  this->registerData(_bottom_loading, &(loads[1]));

  this->registerParameter("shear modulus top", mu[0]);
  this->registerParameter("shear modulus bottom", mu[1]);
  this->registerParameter("poisson ratio top", nu[0]);
  this->registerParameter("poisson ratio bottom", nu[1]);
  this->registerParameter("shear wave speed top", cs[0]);
  this->registerParameter("shear wave speed bottom", cs[1]);
  this->registerParameter("delta x", dx[0]);
  this->registerParameter("delta z", dx[1]);
  this->registerParameter("Domain length X", X[0]);
  this->registerParameter("Domain length Z", X[1]);
  this->registerParameter("beta", beta);
  this->registerParameter("delta min", dxmin);
  
}

/* -------------------------------------------------------------------------- */
void SpectralModel::resetStableTimeStep(Real new_beta) {

  this->beta = new_beta;
  std::cout << "!! Stable time step reinitialized to beta = " << beta << std::endl;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::initConvolutionManagers(bool blank){

  std::stringstream nuk;
  UInt nu_int;
  std::string name1 = "nu_.";
  std::vector<std::string> name3((nb_kernels-1)*2);;
  name3 = {"_h11.dat","_h11b.dat","_h22.dat","_h22b.dat","_k12.dat","_k12b.dat",};

  // Check if the materials are the same on top and bottom, and if so only one occurence of the kernel is required
  
  if((nu[0] == nu[1])&&(mu[0]==mu[1])){
      name3 = {"_h11.dat","_h11.dat","_h22.dat","_h22.dat","_k12.dat","_k12.dat",};
    }

  std::string name2;
 
  std::vector<Real> kside = {1, ksi};
  std::vector<Real> mod_numb;
  Real nu_f;

  SpectralConvolutionManager ** convo_manager;
  Idx size = 0;
 
  for (UInt side = 0; side < 2; ++side) {
    mod_numb.resize(total_nele_fft-1);
    
    for (UInt m = 0; m < total_nele_fft-1; ++m) {
      mod_numb[m] = q00[m]/kside[side];
    }

    if (side==0)
      convo_manager = &convo_manager_top;
    else
      convo_manager = &convo_manager_bot;

    (*convo_manager) = new SpectralConvolutionManager(beta*dxmin, mod_numb, dim, nb_kernels);
    size += (*convo_manager)->init(t_cut[side],ntim,blank);

    for (UInt i = 0; i < (nb_kernels-1); ++i) {
      
      nuk.seekp(std::ios_base::beg);
      nu_int = (UInt)(nu[side]*1000);
      nu_int = 0.1*(nu_int+1);
      nuk << nu_int;
      name2 = nuk.str();
      std::string f90file = name1+name2+name3[2*i+side];
      nu_f = (*convo_manager)->loadKernel(f90file,i);
      
      if (fabs(nu_f - nu[side]) > 1e-4) { 
	std::stringstream err;
	err << "WARNING !!! nu read from " << f90file
	    << " is not consistent with the input parameters"
	    << std::endl;
	cRacklet::error(err);
      }
      std::cout << "Convolution kernel correctly loaded from " << f90file << std::endl;
    }
    ConvolutionManager::KernelFunctor ** funct = (*convo_manager)->getKernelFunctor(3);
    (*funct) = new ConvolutionManager::H_33();
  }
  
  if(blank)
    std::cout << "Required memory size would be: " << size*sizeof(Real)*1e-6 << " MB " << std::endl;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::restartModel(bool from_2dto3d) {

  std::vector<UInt> nfft;
  UInt nx;
  
  if(from_2dto3d) {
    nfft = nele_fft;
    nx = n_ele[0];
  } else {
    nfft = {0,0};
    nx = 0;
  } 
  
  this->restart(false,nx);

  if(interface_law)
    interface_law->restart(false,nx);
  
  UInt step_top = convo_manager_top->restart(0,false,nfft);
  UInt step_bot = convo_manager_bot->restart(1,false,nfft);

  if(step_top==step_bot)
    it = step_top+1;
  else
    cRacklet::error("Mismatching restarting time step between top and bottom material");

  displ_jump->computeJumpFields(); 
  veloc_jump->computeJumpFields();
}

/* -------------------------------------------------------------------------- */
void SpectralModel::pauseModel() {

  std::stringstream dir;
  dir << this->output_dir << restart_dir;
  
  mkdir((dir.str()).c_str(),0777);

  this->restart(true);

  if(interface_law)
    interface_law->restart(true);
  
  convo_manager_top->restart(0,true);
  convo_manager_bot->restart(1,true);
}
/* -------------------------------------------------------------------------- */
void SpectralModel::setLoadingCase(Real load, Real psi, Real phi, bool write) {

  if(write==true)
    printSelfLoad(load, psi, phi);
  
  uniform_loading.resize(dim);
  
  psi *= M_PI/180;
  phi *= M_PI/180;

  uniform_loading[0] = load*sin(psi)*cos(phi);
  
  uniform_loading[1] = load*cos(psi);

  uniform_loading[2] = load*sin(psi)*sin(phi);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::incrementLoad(Real increment, UInt loading_direction) {
  
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      for (UInt side = 0; side < 2; ++side) {
	loads[side][(x*dim+loading_direction)+z*n_ele[0]*dim] += increment;
      }
    }
  }
}


/* -------------------------------------------------------------------------- */
void SpectralModel::updateLoads(Real * loading_per_dim) {
    
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      for (UInt side = 0; side < 2; ++side) {
	for (UInt j = 0; j < dim; ++j) {
	  loads[side][(x*dim+j)+z*n_ele[0]*dim] = *(loading_per_dim+j)*loading_ratio[x+z*n_ele[0]];
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateLoads() {

  updateLoads(&uniform_loading[0]);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::sinusoidalLoading(Real min) {
  
  Real ratio;

  for (UInt x = 0; x < n_ele[0]; ++x) {
    ratio = sin(Real(x)/n_ele[0]*M_PI)*(1-min) + min;
    for (UInt z = 0; z < n_ele[1]; ++z) {
      loading_ratio[x+z*n_ele[0]] = ratio;
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::readSpatialLoadingFromFile(std::string loading_file) {

  std::ifstream fin(loading_file);
  
  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      fin >> loading_ratio[x+z*n_ele[0]];
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::initInterfaceFields() {
  
  veloc_jump->computeJumpFields();
  interface_law->initInterfaceConditions();
  veloc_jump->computeJumpFields();
}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateDisplacements() {

  Real dt = dxmin*beta/(std::max(cs[0],cs[1]));
  
  for (UInt i = 0; i < 2; ++i) {
    displacements[i] += velocities[i]*dt;
  } 
  displ_jump->computeJumpFields(); 
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeInterfaceFields(){

  interface_law->updateInterfaceConditions();
  veloc_jump->computeJumpFields();
}

/* -------------------------------------------------------------------------- */
void SpectralModel::fftOnDisplacements() {

  Real *it_U;
  SpectralConvolutionManager * convo_manager;
  std::vector<Real> U_t(2*(total_nele_fft-1)*dim);
  for (UInt side = 0; side < 2; ++side){
    if (side==0) {
      it_U = U_top;
      convo_manager = convo_manager_top;}
    else {
      it_U = U_bot;
      convo_manager = convo_manager_bot;}
    displacements[side].stridedFFT(it_U, dim);
    convo_manager->storeFields(it_U);
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeStress() {
  
  Real normalize = 1/(Real)total_n_ele; //manuel FFTW

  for (UInt side = 0; side < 2; ++side) {
   
    Real * F_it(F_k);

    if (interface_dim==2)
      computeConvolutions<2>(F_it, side);
    else if (interface_dim==1)
      computeConvolutions<1>(F_it, side);

    for (UInt d = 0; d < dim; ++d) {
      for (UInt img = 0; img < 2; ++img) {
	*(F_k+d*(total_nele_fft*2)+img)=0.0;
      }
    }

    stresses[side].backwardFFT(F_k, dim); 
    stresses[side] *= normalize;
    stresses[side] += loads[side];
    
  }  
}

/* -------------------------------------------------------------------------- */
void SpectralModel::increaseTimeStep() {

 displ_jump->computeJumpFields();
 veloc_jump->computeJumpFields();
 computeAll(it*beta*dxmin/(std::max(cs[0],cs[1])));
 ++it;

 if((it>ntim+1)&&(ntim!=0))
   std::cout << "!!! WARNING! The number of simulation time steps exceed user prediction !!!"
	     << std::endl;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::printSelfLoad(Real load, Real psi, Real phi) {

  out_summary << "/* -------------------------------------------------------------------------- */ " 
	      << std::endl
	      << " LOADING CONDITIONS" << std::endl
	      << "* Loading: " << load << std::endl
	      << "* Psi angle: " << psi << std::endl			      
	      << "* Phi angle: " << phi << std::endl
	      << std::endl;

  out_parameters << "tau_0 " << load << std::endl
		 << "psi " << psi << std::endl
		 << "phi " << phi << std::endl;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::printSelf() {
  
  Real X1;
  if(interface_dim==1){X1=0;}
  else{X1=X[1];}

  out_summary << "/* -------------------------------------------------------------------------- */ "
  << std::endl
  << " MODEL VARIABLES " << std::endl
  << "* Number of elements x: " << n_ele[0] << std::endl
  << "* Number of elements z: " << n_ele[1] << std::endl
  << "* Number of time steps: " << ntim << std::endl
  << "* Stable time step parameter (beta) " << beta << std::endl
  << "* Domain size X: " << X[0] << std::endl
  << "* Domain size Z: " << X1 << std::endl
  << "* Young modulus of top material: " << 2*mu[0]*(1+nu[0]) << std::endl
  << "* Young modulus of bottom material: " << 2*mu[1]*(1+nu[1]) << std::endl
  << "* Poisson ratio of top material: " << nu[0] << std::endl
  << "* Poisson ratio of bottom material: " << nu[1] << std::endl
  << "* Shear wave speed ratio: " << ksi << std::endl
  << std::endl;	

  out_parameters << "nb_ele_x " << n_ele[0] << std::endl
		 << "nb_ele_z " << n_ele[1] << std::endl
		 << "nb_t_step "<< ntim << std::endl
		 << "beta " << beta << std::endl
		 << "dom_size_x " << X[0] << std::endl
		 << "dom_size_z " << X1 << std::endl
		 << "E_top " << 2*mu[0]*(1+nu[0]) << std::endl
		 << "E_bot " << 2*mu[1]*(1+nu[1]) << std::endl
		 << "nu_top "<< nu[0] << std::endl
		 << "nu_bot "<< nu[1] << std::endl
		 << "cs_top " << cs[0] << std::endl
		 << "cs_bot " << cs[1] << std::endl;
}
