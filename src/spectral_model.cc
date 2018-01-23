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
void SpectralModel::initModel(Real reset_beta) {

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

  it = 1;

  nb_kernels = 4;

  zeta = mu[0]/mu[1];

  displacements.resize(2);
  velocities.resize(2);
  stresses.resize(2);
  loads.resize(2);
  eta.resize(2);
 
  for (UInt i = 0; i < 2; ++i) {
    
    displacements[i].SetGridSize(n_ele, dim);
    velocities[i].SetGridSize(n_ele, dim);
    stresses[i].SetGridSize(n_ele, dim);
    loads[i].SetGridSize(n_ele, dim);
    eta[i] = sqrt(2*(1-nu[i])/(1-2*nu[i]));

    displacements[i].initFFT(true,dim);
    stresses[i].initFFT(false,dim);
  }

  displ_jump = new InterfaceFields(&displacements, dim);
  veloc_jump = new InterfaceFields(&velocities, dim);

  intfc_trac.SetGridSize(n_ele, dim);
  nor_opening.SetGridSize(n_ele, 1);
  shr_opening.SetGridSize(n_ele, 1);

  nor_strength.resize(total_n_ele);
  shr_strength.resize(total_n_ele);
  fric_strength.resize(total_n_ele);
  ind_crack.resize(total_n_ele);
 
  /* -------------------------------------------------------------------------- */
  //set the modes q0, k and m
  if(interface_dim==2)
    initFrequency<2>();
  else if(interface_dim==1)
    initFrequency<1>();

  h11u1 = new Real[2*(total_nele_fft-1)];
  h22u2 = new Real[2*(total_nele_fft-1)];
  h33u3 = new Real[2*(total_nele_fft-1)];
  h12u1 = new Real[2*(total_nele_fft-1)];
  h12u2 = new Real[2*(total_nele_fft-1)];

  U_top = new Real[2*(total_nele_fft-1)*dim];
  U_bot = new Real[2*(total_nele_fft-1)*dim];
  F_k = new Real[2*total_nele_fft*dim];
  
  this->nb_kernels = 4;

  initConvolutionManagers();

  registerModelFields();
  
  printSelf();
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
  this->registerData(_interface_tractions, &intfc_trac);
  this->registerData(_top_loading, &(loads[0]));
  this->registerData(_bottom_loading, &(loads[1]));
  this->registerData(_normal_strength, &nor_strength);
  this->registerData(_shear_strength, &shr_strength);
  this->registerData(_frictional_strength, &fric_strength);
  this->registerData(_id_crack, &ind_crack);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::resetStableTimeStep(Real new_beta) {

  this->beta = new_beta;
  std::cout << "!! Stable time step reinitialized to beta = " << beta << std::endl;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::initConvolutionManagers(){

  std::stringstream nuk;
  UInt nu_int;
  std::string name1 = "nu_.";
  std::vector<std::string> name3((nb_kernels-1)*2);;
  name3 = {"_h11.dat","_h11b.dat","_h22.dat","_h22b.dat","_k12.dat","_k12b.dat",};
  std::string name2;
 
  std::vector<Real> kside = {1, ksi};
  std::vector<Real> mod_numb;
  Real nu_f;

  SpectralConvolutionManager ** convo_manager;

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
    (*convo_manager)->init(t_cut[side],ntim);

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

  if(contact_law)
    contact_law->restart(false,nx);
  if(fracture_law)
    fracture_law->restart(false,nx);
  
  UInt step_top = convo_manager_top->restart(0,false,nfft);
  UInt step_bot = convo_manager_bot->restart(1,false,nfft);

  if(step_top==step_bot)
    it = step_top+1;
  else
    cRacklet::error("Mismatching restarting time step between top and bottom material");
}

/* -------------------------------------------------------------------------- */
void SpectralModel::pauseModel() {

  std::stringstream dir;
  dir << this->output_dir << restart_dir;
  
  mkdir((dir.str()).c_str(),0777);

  this->restart(true);

  if(contact_law)
    contact_law->restart(true);
  if(fracture_law)
    fracture_law->restart(true);
  
  convo_manager_top->restart(0,true);
  convo_manager_bot->restart(1,true);
}
/* -------------------------------------------------------------------------- */
void SpectralModel::setLoadingCase(Real load, Real psi, Real phi) {

  printSelfLoad(load, psi, phi);
  
  uniform_loading.resize(dim);
  
  psi *= M_PI/180;
  phi *= M_PI/180;

  uniform_loading[0] = load*sin(psi)*cos(phi);
  
  uniform_loading[1] = load*cos(psi);

  uniform_loading[2] = load*sin(psi)*sin(phi);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateLoads(Real * loading_per_dim) {

  for (UInt x = 0; x < n_ele[0]; ++x) {
    for (UInt z = 0; z < n_ele[1]; ++z) {
      for (UInt side = 0; side < 2; ++side) {
	for (UInt j = 0; j < dim; ++j) {
	  loads[side][(x*dim+j)+z*n_ele[0]*dim] = *(loading_per_dim+j);     
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
  
  loading_ratio = new Real[total_n_ele];
  Real ratio;

  for (UInt x = 0; x < n_ele[0]; ++x) {
    ratio = sin(Real(x)/n_ele[0]*M_PI)*(1-min) + min;
    for (UInt z = 0; z < n_ele[1]; ++z) {
      loading_ratio[x+z*n_ele[0]] = ratio;
    }
  }
}

/* -------------------------------------------------------------------------- */
#ifdef CRACKLET_USE_LIBSURFER 
void SpectralModel::brownianHeterogLoading(Real rms, long int seed,
					   Real hurst, UInt q0,
					   UInt q1, UInt q2){
  
  SurfaceGeneratorFilterFFT surf_gen;
  int & grid_size = surf_gen.getGridSize();
  grid_size = std::max(4*n_ele[0],4*n_ele[1]);
  Real & Hurst = surf_gen.getHurst();
  Hurst = hurst;
  Real & RMS = surf_gen.getRMS();
  RMS = rms;
  int & Q0 = surf_gen.getQ0();
  Q0 = q0;
  int & Q1 = surf_gen.getQ1();
  Q1 = q1;
  int & Q2 = surf_gen.getQ2();
  Q2 = q2;
  long int & Seed = surf_gen.getRandomSeed();
  Seed = seed;
  surf_gen.Init();
  Surface<Real> & surface = surf_gen.buildSurface();
  std::cout << "Successfully generated loading distribution with an RMS of " << SurfaceStatistics::computeStdev(surface) << std::endl;

  loading_ratio = new Real[total_n_ele];

  for (UInt ix = 0; ix < n_ele[0]; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {
      loading_ratio[ix+iz*n_ele[0]] = surface(ix,iz);      
    }
  }
}
#endif

/* -------------------------------------------------------------------------- */
void SpectralModel::computeInitialVelocities() {

  Real strength;
  Real shr_trac;
  Real shr_velo;
  std::vector<Real> temp_f(2);
  
  UInt i=0;
  for (UInt h = 0; h < n_ele[0]; ++h) {
    for (UInt j = 0; j < n_ele[1]; ++j) {
      i=h+j*n_ele[0];
    
      if((nor_strength[i]==0)&&(loads[0][i*dim+1] < 0.0)) { 
      
	contact_law->computeFricStrength(loads[0][i*dim+1], strength, i, it); 
      
	for (UInt side = 0; side < 2; ++side) {
	  velocities[side][i*dim+1] = 0.0;
	}
      }
      else{//velocities u2
	strength = shr_strength[i];
	velocities[0][i*dim+1] = std::max((loads[0][i*dim+1]-nor_strength[i])/(mu[0]*eta[0]),0.0);
	velocities[1][i*dim+1] = std::min((zeta/ksi)*(nor_strength[i]-loads[1][i*dim+1])/(mu[0]*eta[1]),0.0);
      }

      //velocities u1 & u3
      for (UInt side = 0; side < 2; ++side) {
  
	for (UInt k = 0; k < 2; ++k) {
	  temp_f[k] = loads[side][i*dim+2*k];
	}
      
	shr_trac = sqrt(temp_f[0]*temp_f[0]+temp_f[1]*temp_f[1]);
	if (shr_trac ==0) {
	
	  for (UInt k = 0; k < 2; ++k) {
	  
	    velocities[side][i*dim+2*k] = 0.0;
	  }
	}
	else {
	  if (side==0) shr_velo = std::max((shr_trac-strength)/mu[0],0.0);
	  else shr_velo = std::min((zeta/ksi)*(strength - shr_trac)/mu[0],0.0);
	  
	  for (UInt k = 0; k < 2; ++k) {
	    velocities[side][i*dim+2*k] = shr_velo*temp_f[k]/shr_trac;
	  } 
	}
      }
    }
  }
  veloc_jump->computeJumpFields();
}
/* -------------------------------------------------------------------------- */
void SpectralModel::updateDisplacements() {

  for (UInt i = 0; i < 2; ++i) {
    displacements[i] += velocities[i]*dxmin*beta;
  }

  displ_jump->computeJumpFields(); 

  shr_opening = displ_jump->delta_fields[0];
  nor_opening = displ_jump->delta_fields[1];

  }

/* -------------------------------------------------------------------------- */
void SpectralModel::updateMaterialProp(){

 
  fracture_law->updateFractureLaw(nor_strength, shr_strength, ind_crack, nor_opening, shr_opening);  
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
 void SpectralModel::computeVelocities(){

   CrackProfile deltaStresses(n_ele,dim);
   std::vector<Real> temp_veloc(dim);
   Real trac;
   
   deltaStresses = stresses[0] - stresses[1];
   
   Real cste = 1/(mu[0]*(1+ksi/zeta)); 

   velocities[0] =  deltaStresses * cste;

   for (UInt i = 0; i < n_ele[0]; ++i) {
     for (UInt j = 0; j < n_ele[1]; ++j) {
        velocities[0][(i*dim+1)+j*n_ele[0]*dim] *= (1+ksi/zeta)/(eta[0]+ksi*eta[1]/zeta); 
     }
   }
   velocities[1]=velocities[0]; 

   for (UInt i = 0; i < n_ele[0]; ++i) {
     for (UInt j = 0; j < n_ele[1]; ++j) {

       trac = stresses[0][(i*dim+1)+j*n_ele[0]*dim] - mu[0]*eta[0]*velocities[0][(i*dim+1)+j*n_ele[0]*dim];
       if ((nor_strength[i+n_ele[0]*j] < trac)||(nor_strength[i+n_ele[0]*j]==0)) computeIndepNormalVelocities(i,j);
       else {
    
	 intfc_trac[(i*dim+1)+j*n_ele[0]*dim] = trac;
	 computeShearVelocities(shr_strength[i + n_ele[0]*j], i + j*n_ele[0]);
       }
     }
   }
   veloc_jump->computeJumpFields();
 }

/* -------------------------------------------------------------------------- */
void SpectralModel::computeIndepNormalVelocities(UInt ix, UInt iz){

  std::vector<Real> temp_veloc(2);
  std::vector<Real> cmpted_stress(2);
  Real delta_overlap;

  UInt i = ix+iz*n_ele[0];

  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = stresses[side][i*dim+1]; 
  }

  temp_veloc[0] = 1/(mu[0]*eta[0]) * (cmpted_stress[0]-nor_strength[i]); 
  temp_veloc[1] = 1/(mu[0]*eta[1]) * (zeta/ksi) * (nor_strength[i] - cmpted_stress[1]);

  delta_overlap = displacements[0][i*dim+1] - displacements[1][i*dim+1] + beta*dxmin*(temp_veloc[0]-temp_veloc[1]); 

  if (cRacklet::is_negative(delta_overlap)&&(!overlapping)) computeContactVelocities(ix, iz); 
  else {

    for (UInt side = 0; side < 2; ++side) {

      velocities[side][i*dim+1] = temp_veloc[side];
    }
    intfc_trac[i*dim+1] = nor_strength[i]; 
    computeShearVelocities(shr_strength[i], i);
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeContactVelocities(UInt ix, UInt iz){

  Real aux;
  Real temp_velot;
  Real temp_trac;
  Real strength;
  std::vector<Real> cmpted_stress(2);

  UInt i = ix+iz*n_ele[0];

  for (UInt side = 0; side < 2; ++side) {
    cmpted_stress[side] = stresses[side][i*dim+1];
  }

  aux = (displacements[0][i*dim+1] - displacements[1][i*dim+1])/(dxmin*beta);
  temp_velot = 1/(eta[0]+ksi*eta[1]/zeta)*((cmpted_stress[0]-cmpted_stress[1])/mu[0] - aux*ksi*eta[1]/zeta);
  velocities[0][i*dim+1] = temp_velot;
  velocities[1][i*dim+1] = temp_velot + aux;

  temp_trac = cmpted_stress[0] - eta[0]*mu[0]*temp_velot;
  intfc_trac[i*dim+1] = temp_trac;

  contact_law->computeFricStrength(temp_trac, strength, i, it);

  computeShearVelocities(strength, i);
  
  ind_crack[i] = 3;
  fric_strength[i] = strength;
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeShearVelocities(Real strength, UInt i) {

  std::vector<Real> trac(2);
  Real shr_trac;

  for (UInt j = 0; j < 2; ++j) {
  
    trac[j] = stresses[0][i*dim+2*j] - mu[0]*velocities[0][i*dim+2*j]; 
  }
 
  shr_trac = sqrt((trac[0]*trac[0])+(trac[1]*trac[1])); 

  if ((strength < shr_trac)||(strength==0)) computeIndepShearVelocities(strength, i);
  else{
    
    for (UInt j = 0; j < (dim-1); ++j) {
      
      intfc_trac[i*dim+2*j] = trac[j];
    }
  }

}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeIndepShearVelocities(Real strength, UInt i){

   std::vector<Real> cmpted_stress(2);
   Real dyn_stress;
   Real shr_veloc;

   for (UInt side = 0; side < 2; ++side) {
     
     for (UInt j = 0; j < 2; ++j) {
       
       cmpted_stress[j] = stresses[side][i*dim+2*j];
     }

     dyn_stress = sqrt((cmpted_stress[0]*cmpted_stress[0])+(cmpted_stress[1]*cmpted_stress[1])); 

     if (side==0) shr_veloc = 1/mu[0]*(dyn_stress-strength); 
     
     else shr_veloc = 1/mu[0]*(zeta/ksi)*(strength-dyn_stress); 

      

     for (UInt j = 0; j < 2; ++j) {
       
       if(dyn_stress==0){velocities[side][i*dim+2*j]=0;}
       else{velocities[side][i*dim+2*j] = shr_veloc*cmpted_stress[j]/dyn_stress;}
       if (side==0){
	 if(dyn_stress==0){intfc_trac[i*dim+2*j] =0;}
	 else{intfc_trac[i*dim+2*j] = strength*cmpted_stress[j]/dyn_stress;}
       }
     }    
   }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::increaseTimeStep() {

 displ_jump->computeJumpFields();
 veloc_jump->computeJumpFields();
 computeAll(it*beta*dxmin/X[0]);
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
		 << "cs_top " << cs_t << std::endl
		 << "cs_bot " << cs_t/ksi << std::endl;
}
