/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include <vector>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdint>
#include <gsl/gsl_sf_bessel.h>
#include <complex>
#include <fftw3.h>
#include <math.h>
#include <algorithm>

#define PI 3.141592653589

/* -------------------------------------------------------------------------- */
void Energetics::integrate(double dt){

  E += 0.5*dt*(E_dot_old + E_dot);
  E_dot_old = E_dot;
  E_dot = 0;

}

/* -------------------------------------------------------------------------- */
void InterfaceFields::computeJumpFields() {

  delta_fields.resize(2);

  CrackProfile delta_field_1 = (*fields)[0].getStridedPart(dim,0) - (*fields)[1].getStridedPart(dim,0);
  CrackProfile delta_field_3 = (*fields)[0].getStridedPart(dim,2) - (*fields)[1].getStridedPart(dim,2);

  delta_fields[0] = delta_field_1*delta_field_1 + delta_field_3*delta_field_3;
  delta_fields[0].squareRoot();

  delta_fields[1] = (*fields)[0].getStridedPart(dim,1) - (*fields)[1].getStridedPart(dim,1);
}

/* -------------------------------------------------------------------------- */
void SpectralModel::initModel() {

  dim = 3;

  dx = X/(double)(n_ele);

  q0 = 2*PI/(dx*(double)n_ele);

  it = 1;

  nele_fft = n_ele/2;

  beta = 0.4;

  nb_kernels = 4;

  zeta = mu[0]/mu[1];

  displacements.resize(2);
  velocities.resize(2);
  stresses.resize(2);
  loads.resize(2);
  kernel.resize(2);
  eta.resize(2);
  t_cut_j.resize(2);

  E_nor.resize(2);
  E_shr.resize(2);
  E_fri.resize(2);

  for (int i = 0; i < 2; ++i) {
    
    displacements[i].SetGridSize(n_ele*dim);
    velocities[i].SetGridSize(n_ele*dim);
    stresses[i].SetGridSize(n_ele*dim);
    loads[i].SetGridSize(n_ele*dim);
    eta[i] = sqrt(2*(1-nu[i])/(1-2*nu[i]));

  }

  displ_jump = new InterfaceFields(&displacements, dim);
  veloc_jump = new InterfaceFields(&velocities, dim);

  intfc_trac.SetGridSize(n_ele*dim);
  nor_opening.SetGridSize(n_ele);
  shr_opening.SetGridSize(n_ele);

  nor_strength.resize(n_ele);
  shr_strength.resize(n_ele);
  ind_crack.resize(n_ele);

  U.resize(2*dim);
  K.resize(2*nb_kernels); 
  convo.resize(2*5);

  for (int i = 0; i < 2*nb_kernels; ++i) {
    K[i].resize(nele_fft);
  }

for (int i = 0; i < 2*5; ++i) {
  convo[i].resize(nele_fft);
 }


  t_cut_j.resize(2);
  
  for (int side = 0; side < 2; ++side) {

    t_cut_j[side].resize(nele_fft);

    for (int j = 1; j <= nele_fft; ++j) {

      if(side==0) t_cut_j[side][(j-1)] = std::min(ntim,
						  (int)(t_cut[0]*(double)n_ele*0.5/(PI*(double)j*beta))-1);
      if(side==1) t_cut_j[side][(j-1)] = std::min(ntim,
						  (int)(ksi*t_cut[1]*(double)n_ele*0.5/(PI*(double)j*beta))-1);
   
      for (int i = 0; i < nb_kernels; ++i) {
	K[2*i+side][j-1].resize(t_cut_j[side][(j-1)]);
      }
    }
    
    for (int i = 0; i < dim; ++i) {
      U[2*i+side].resize(t_cut_j[side]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::setLoadingCase(double load_in, double load_out, double psi, double phi) {

  loads_insd.resize(3);
  loads_outsd.resize(3);

  psi *= PI/180;
  phi *= PI/180;

  loads_insd[0] = load_in*sin(psi)*cos(phi);
  loads_outsd[0] = load_out*sin(psi)*cos(phi);

  loads_insd[1] = load_in*cos(psi);
  loads_outsd[1] = load_out*cos(psi);

  loads_insd[2] = load_in*sin(psi)*sin(phi);
  loads_outsd[2] = load_out*sin(psi)*sin(phi);

}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateLoads() {


  for (int i = 0; i < n_ele; ++i) {
     
    if (ind_crack[i]>=load_index) {
      for (int side = 0; side < 2; ++side) {
	for (int j = 0; j < dim; ++j) {
	  loads[side][i*dim+j] = loads_insd[j];     
	}
      }
    }
    else {
      for (int side = 0; side < 2; ++side) {
	for (int j = 0; j < dim; ++j) {
	  loads[side][i*dim+j] = loads_outsd[j];
	}
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateTriggerLoading(double alpha, double beta, double ampl) {

  for (int i = 0; i < n_ele; ++i) {

    for (int side = 0; side < 2; ++side) {

      double delta = (1+beta)*ampl/cosh(alpha*((double)i/n_ele - 0.5));

      if (delta > loads_insd[0]){

	loads[side][i*dim] = delta;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeInitialVelocities() {

  double strength;
  double shr_trac;
  double shr_velo;
  double dv1, dv3, dvn, dvs;
  std::vector<double> temp_f(2);
  
  for (int i = 0; i < n_ele; ++i) {
    
    if((nor_strength[i]==0)&&(loads[0][i*dim+1] < 0)) {
      
      contact_law->computeFricStrength(loads[0][i*dim+1], strength, i, it);
      
      for (int side = 0; side < 2; ++side) {
	velocities[side][i*dim+1] = 0.0;
      }
    }

    else{

      strength = shr_strength[i];
      velocities[0][i*dim+1] = std::max((loads[0][i*dim+1]-nor_strength[i])/(mu[0]*eta[0]),0.0);
      velocities[1][i*dim+1] = std::min((zeta/ksi)*(nor_strength[i]-loads[1][i*dim+1])/(mu[0]*eta[1]),0.0);
    }
    
    for (int side = 0; side < 2; ++side) {
  
      for (int k = 0; k < 2; ++k) {
	temp_f[k] = loads[side][i*dim+2*k];
      }
      
      shr_trac = sqrt(temp_f[0]*temp_f[0]+temp_f[1]*temp_f[1]);
      
      if (shr_trac ==0) {
	
	for (int k = 0; k < 2; ++k) {
	  
	  velocities[side][i*dim+2*k] = 0.0;
	}
      }
      else {
	
	if (side==0) shr_velo = std::max((shr_trac-strength)/mu[0],0.0);
	else shr_velo = std::min((zeta/ksi)*(strength - shr_trac)/mu[0],0.0);
	  
	for (int k = 0; k < 2; ++k) {
	  velocities[side][i*dim+2*k] = shr_velo*temp_f[k]/shr_trac;
	} 
      }
    }

    dv1 = velocities[0][i*dim] - velocities[1][i*dim];
    dv3 = velocities[0][i*dim+2] - velocities[1][i*dim+2];
    dvs = sqrt(dv1*dv1 + dv3*dv3);
    dvn = velocities[0][i*dim+1] - velocities[1][i*dim+1];

    if (i<n_ele/2) {
      E_nor[0].E_dot +=  nor_strength[i]*fabs(dvn)*dx;
      E_shr[0].E_dot += shr_strength[i]*fabs(dvs)*dx;
    }
    else {
      E_nor[1].E_dot +=  nor_strength[i]*fabs(dvn)*dx;
      E_shr[1].E_dot += shr_strength[i]*fabs(dvs)*dx;
    }

    if((nor_strength[i]==0)&&(loads[0][i*dim+1] < 0)) {
      
      if (i<n_ele/2)	 E_fri[0].E_dot += strength * fabs(dvs) * dx;
      else             E_fri[1].E_dot += strength * fabs(dvs) * dx;
    }
  }

  for (int i = 0; i < 2; ++i) {
    
    E_nor[i].integrate(beta*dx);
    E_shr[i].integrate(beta*dx);
    E_fri[i].integrate(beta*dx);
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateDisplacements() {

  for (int i = 0; i < 2; ++i) {
    displacements[i] = displacements[i]+ velocities[i]*dx*beta;
  }

  displ_jump->computeJumpFields();

  shr_opening = (*displ_jump).delta_fields[0];
  nor_opening = (*displ_jump).delta_fields[1];

}

/* -------------------------------------------------------------------------- */
void SpectralModel::updateMaterialProp(){

 
  fracture_law->updateFractureLaw(nor_strength, shr_strength, ind_crack, nor_opening, shr_opening);

  updateLoads();

}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeKernels(){

  bool side_index;
  
  delta_smpl.resize((nb_kernels-1)*2);
  kernel.resize((nb_kernels-1)*2);

  std::stringstream nuk;
  int nu_int;
  std::string name1 = "nu_.";
  std::vector<std::string> name3((nb_kernels-1)*2);;
  name3 = {"_h11.dat","_h11b.dat","_k12.dat","_k12b.dat","_h22.dat","_h22b.dat"};
  std::string name2;
   
  for (int i = 0; i < (nb_kernels-1); ++i) {

    for (int side = 0; side < 2; ++side) {
      
    nuk.seekp(std::ios_base::beg);
    nu_int = (int)(nu[side]*1000);
    nu_int = 0.1*(nu_int+1);
    nuk << nu_int;
    name2 = nuk.str();
    readKernel(name1+name2+name3[2*i+side],2*i+side);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::readKernel(std::string f90file, int ik){


  std::ifstream file (f90file, std::ios::binary);

  if (!file.is_open()){
    std::cout << "Unable to open file" << std::endl;
    std::cout << "Check that the file " << f90file 
	      <<" is in the current folder" << std::endl;
    return;
  }
 
  int inu;
  int side_index;

  side_index = ((double)(ik/2)==(double)ik/2);

  if (side_index) inu=0;
  else inu=1;
  
  double temp;
  double cdcs_f;
  double nu_f;
  int numt;

  readBinary(file, numt);
  readBinary(file, temp);
  if (ik==0) t_cut[0] = temp;
  else if (ik == 1) t_cut[1] = temp;
  readBinary(file, numt);
  readBinary(file, delta_smpl[ik]);
  readBinary(file, cdcs_f);
  readBinary(file, nu_f);
 
  kernel[ik].resize(numt);
  readBinary(file, kernel[ik]);
  
  // Test

  if (fabs(cdcs_f - eta[inu]) > 1e-4) std::cout << "WARNING ! cd/cs read from " << f90file
						<< " is not consistent with the input parameters" 
						<< std::endl;
  if (fabs(nu_f - nu[inu]) > 1e-4) std::cout << "WARNING ! nu read from " << f90file
					    << " is not consistent with the input parameters"
					    << std::endl;

  std::cout << "Convolution kernel correctly loaded from " << f90file << std::endl;

  for (int i = 0; i < numt; ++i) {

    if (std::isnan(kernel[ik][i])) {

      std::cout << "!! NaN values read for kernel n° " << ik << ", " << i << "th values !!" << std::endl;
    
  }
  
 }

}

/* -------------------------------------------------------------------------- */
double SpectralModel::H33(double x) {

  double value;

  value = gsl_sf_bessel_J1(x)/x;

  return value;

}

/* -------------------------------------------------------------------------- */
void SpectralModel::preintegratedKernels() {

  double up_bound0;
  double low_bound0;

  double up_bound;
  double low_bound;

  double trap_sum;
  double nu_;
  bool side_index;
  bool cut_top;
  bool cut_bot;

  up_bound0 = q0*dx*beta* (double)it;
  low_bound0 = q0*dx*beta* ((double)it-1+1e-10);

  for (int j = 1; j <= nele_fft; ++j) {

    cut_top =  it > t_cut_j[0][j-1];
    cut_bot = it > t_cut_j[1][j-1];
    
    if ((cut_top) && (cut_bot)) break;
    
    for (int i = 0; i < nb_kernels; ++i) {
     
      for (int side = 0; side < 2; ++side) {
           
	if ((side==0) && (!cut_top)) {
	  up_bound = up_bound0 * (double)j;
	  low_bound = low_bound0 * (double)j;
	}
	else if ((side==1) && (!cut_bot)) {
	  up_bound = up_bound0 * (double)j / ksi;
	  low_bound = low_bound0 * (double)j / ksi;
	}
	else continue;
     
	if (i < (nb_kernels-1)) { 

	  trap_sum = 0.5 * (interpolateFromKernel(up_bound,2*i+side) 
			+ interpolateFromKernel(low_bound,2*i+side));
	}
	else {
	  trap_sum = 0.5 * (H33(up_bound) + H33(low_bound));
	}

	K[2*i+side][j-1][it-1] = mu[side] * q0 * (double)j * (up_bound-low_bound) * trap_sum;

      }
    }
  }
}

/* -------------------------------------------------------------------------- */
double SpectralModel::interpolateFromKernel(double val, int index) {

  int i_low;
  double val_low;
  double val_up;
  double ker_low;
  double ker_up;
  double result;
  double delta = delta_smpl[index];

  i_low = (int)(val/delta);

  ker_low = kernel[index][i_low];
  ker_up = kernel[index][i_low+1];

  val_low = ((double) i_low) * delta;
  val_up = ((double) i_low +1)*delta;

  result = (ker_up-ker_low) * (val-val_low)/(val_up-val_low) + ker_low;

  return result;
  
}

/* -------------------------------------------------------------------------- */
void SpectralModel::fftOnDisplacements() {

  std::vector<double> U_t(2*nele_fft*dim);

  for (int side = 0; side < 2; ++side){

    displacements[side].stridedFFT(U_t, dim);

    for (int i = 0; i < dim; ++i) {
      
      U[2*i+side].store(&U_t[i*nele_fft*2] , it-1);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeTimeConvolution() {


  for (int side = 0; side < 2; ++side) {
    
    for (int i = 0; i < dim; ++i) {  
	
      if (i < 2) U[2*i+side].convolution(K[2*i+side], convo[2*i+side], it);
	
      U[2*i+side].convolution(K[2*(i+1)+side], convo[2*i+side+4], it); 
    }
  } 
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeStresses() {

  std::vector<std::complex<double> > U_1(nele_fft);
  std::vector<std::complex<double> > U_2(nele_fft);
  std::complex<double> i;
  double normalize = 1/(double)n_ele;
  std::vector<double> sgn = {1.0, -1.0};
  std::vector<std::complex<double> > F(dim*nele_fft);
 
  i = {0.0, 1.0};

  for (int side = 0; side < 2; ++side) {

    F.resize(dim*(nele_fft+1));

    for (int d = 0; d < dim; ++d) {

      F[d*nele_fft]={0.0 , 0.0}; 

    }
  
    U[side].recall(U_1, it-1);
    U[2+side].recall(U_2, it-1);

    for (int j = 0; j < nele_fft; ++j) {
      
      F[j+1] = -sgn[side] * convo[side][j] + i * convo[2+side][j] 
	+ i*(2-eta[side])*mu[side]*q0*(double)(j+1)*U_2[j];
      
      F[j+nele_fft+2] = -sgn[side] * convo[6+side][j] - i * convo[4+side][j] 
	- i*(2-eta[side])*mu[side]*q0*(double)(j+1)*U_1[j];
      
      F[j+2*nele_fft+3] = -sgn[side] * convo[8+side][j]; 
      
    }

    stresses[side].backwardFFT(&F[0], dim); 
   stresses[side] = loads[side] + stresses[side] * normalize;

  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeVelocities(){

  CrackProfile deltaStresses(n_ele*dim);
  std::vector<double> temp_veloc(dim);
  double trac;

  deltaStresses = stresses[0] - stresses[1];

  double cste = 1/(mu[0]*(1+ksi/zeta)); 

  velocities[0] =  deltaStresses * cste;

  for (int i = 0; i < n_ele; ++i) {
 
    velocities[0][i*dim+1] *= (1+ksi/zeta)/(eta[0]+ksi*eta[1]/zeta); 
  }
  velocities[1]=velocities[0];

  for (int i = 0; i < n_ele; ++i) {

    trac = stresses[0][i*dim+1] - mu[0]*eta[0]*velocities[0][i*dim+1];
    if ((nor_strength[i] < trac)||(nor_strength[i]==0)) computeIndepNormalVelocities(i);
    else {
    
      intfc_trac[i*dim+1] = trac;
      computeShearVelocities(shr_strength[i], i);
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeIndepNormalVelocities(int i){

  std::vector<double> temp_veloc(2);
  std::vector<double> cmpted_stress(2);
  double delta_overlap;

  for (int side = 0; side < 2; ++side) {
    cmpted_stress[side] = stresses[side][i*dim+1];
  }

  temp_veloc[0] = 1/(mu[0]*eta[0]) * (cmpted_stress[0]-nor_strength[i]);
  temp_veloc[1] = 1/(mu[0]*eta[1]) * (zeta/ksi) * (nor_strength[i] - cmpted_stress[1]);

  delta_overlap = displacements[0][i*dim+1] - displacements[1][i*dim+1] + beta*dx*(temp_veloc[0]-temp_veloc[1]);

  if ((delta_overlap < 0)&&(!overlapping)) computeContactVelocities(i);
  else {

    for (int side = 0; side < 2; ++side) {

      velocities[side][i*dim+1] = temp_veloc[side];
    }
    intfc_trac[i*dim+1] = nor_strength[i]; 
    computeShearVelocities(shr_strength[i], i);
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeContactVelocities(int i){

  double aux;
  double temp_velot;
  double temp_trac;
  double strength;
  std::vector<double> cmpted_stress(2);
  double dvs_top, dvs_bot, dvs;

  for (int side = 0; side < 2; ++side) {
    cmpted_stress[side] = stresses[side][i*dim+1];
  }

  aux = (displacements[0][i*dim+1] - displacements[1][i*dim+1])/(dx*beta);

  temp_velot = 1/(eta[0]+ksi*eta[1]/zeta)*((cmpted_stress[0]-cmpted_stress[1])/mu[0] - aux*ksi*eta[1]/zeta);
  velocities[0][i*dim+1] = temp_velot;
  velocities[1][i*dim+1] = temp_velot + aux;

  temp_trac = cmpted_stress[0] - eta[0]*mu[0]*temp_velot;
  intfc_trac[i*dim+1] = temp_trac;

  contact_law->computeFricStrength(temp_trac, strength, i, it);

  computeShearVelocities(strength, i);
  
  dvs_top = sqrt(velocities[0][i*dim]*velocities[0][i*dim] + velocities[0][i*dim+2]*velocities[0][i*dim+2]);
  dvs_bot = sqrt(velocities[1][i*dim]*velocities[1][i*dim] + velocities[1][i*dim+2]*velocities[1][i*dim+2]);
  dvs = dvs_top - dvs_bot;

  if (i<n_ele/2)  E_fri[0].E_dot += strength * fabs(dvs) * dx;
  else            E_fri[1].E_dot += strength * fabs(dvs) * dx; 

}
/* -------------------------------------------------------------------------- */
void SpectralModel::computeShearVelocities(double strength, int i) {

  std::vector<double> trac(2);
  double shr_trac;

  for (int j = 0; j < 2; ++j) {
  
    trac[j] = stresses[0][i*dim+2*j] - mu[0]*velocities[0][i*dim+2*j];
  }
 
  shr_trac = sqrt((trac[0]*trac[0])+(trac[1]*trac[1])); 

  if ((strength < shr_trac)||(strength==0)) computeIndepShearVelocities(strength, i);
  else{
    
    for (int j = 0; j < (dim-1); ++j) {
      
      intfc_trac[i*dim+2*j] = trac[j];
    }
  }

}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeIndepShearVelocities(double strength, int i){

   std::vector<double> cmpted_stress(2);
   double dyn_stress;
   double shr_veloc;

   for (int side = 0; side < 2; ++side) {
     
     for (int j = 0; j < 2; ++j) {
       
       cmpted_stress[j] = stresses[side][i*dim+2*j];
     }

     dyn_stress = sqrt((cmpted_stress[0]*cmpted_stress[0])+(cmpted_stress[1]*cmpted_stress[1]));

     if (side==0) shr_veloc = 1/mu[0]*(dyn_stress-strength);
     
     else shr_veloc = 1/mu[0]*(zeta/ksi)*(strength-dyn_stress);

      

     for (int j = 0; j < 2; ++j) {
       
       velocities[side][i*dim+2*j] = shr_veloc*cmpted_stress[j]/dyn_stress;
       if (side==0) intfc_trac[i*dim+2*j] = strength*cmpted_stress[j]/dyn_stress;
     }    
   }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::computeEnergy() {

  veloc_jump->computeJumpFields();

  CrackProfile delta_vs = (*veloc_jump).delta_fields[0];

  CrackProfile delta_vn = (*veloc_jump).delta_fields[1];

  for (int i = 0; i < n_ele/2; ++i) {

    E_nor[0].E_dot += nor_strength[i]*fabs(delta_vn[i])*dx;
    E_shr[0].E_dot += shr_strength[i]*fabs(delta_vs[i])*dx;
    
    E_nor[1].E_dot += nor_strength[i+n_ele/2]*fabs(delta_vn[i+n_ele/2])*dx;
    E_shr[1].E_dot += shr_strength[i+n_ele/2]*fabs(delta_vs[i+n_ele/2])*dx;
  }

  for (int i = 0; i < 2; ++i) {
    
    E_nor[i].integrate(beta*dx);
    E_shr[i].integrate(beta*dx);
    E_fri[i].integrate(beta*dx);
    
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::generatePulse(double delta, double lx) {

  int nb_elem = (int)(lx/dx);

  int start = 0.5*(n_ele-nb_elem);
  int stop = 0.5*(n_ele+nb_elem);

  for (int i = start; i < stop; ++i) {
  
    for (int side = 0; side < 2; ++side) {

      loads[side][i*dim] *= delta;
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralModel::printSelf(std::ofstream & parameters_file, std::ofstream & summary) {

  double load_in_magn, load_out_magn, psi, phi;

  load_in_magn = 0;
  load_out_magn = 0;

  for (int i = 0; i < dim; ++i) {

    load_in_magn += loads_insd[i]*loads_insd[i];
    load_out_magn += loads_outsd[i]*loads_outsd[i];
  
  }

  psi = acos(loads_insd[1]/sqrt(load_in_magn));
  if (psi!=0) phi = acos(loads_insd[0]/(sin(psi)*sqrt(load_in_magn)));
  else phi = asin(loads_insd[2]/(sin(psi)*sqrt(load_in_magn)));
  

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " MODEL VARIABLES " << std::endl;
  summary << "* Number of elements: " << n_ele << std::endl;		        
  summary << "* Number of time steps: " << ntim << std::endl;		      
  summary << "* Domain size: " << X << std::endl;			      
  summary << "* Initial crack size: " << a_0 << std::endl;		
  summary << "* Load magnitude (inside crack): " << sqrt(load_in_magn) << std::endl;
  summary << "* Load magnitude (bonded interface): " << sqrt(load_out_magn) << std::endl;
  summary << "* Psi angle: " << 180*psi/PI << std::endl;			      
  summary << "* Phi angle: " << 180*phi/PI << std::endl;
  summary << "* Young modulus of top material: " << 2*mu[0]*(1+nu[0]) << std::endl;
  summary << "* Young modulus of bottom material: " << 2*mu[1]*(1+nu[1]) << std::endl;
  summary << "* Poisson ratio of top material: " << nu[0] << std::endl;	      
  summary << "* Poisson ratio of bottom material: " << nu[1] << std::endl;    
  summary << "* Shear wave speed ratio: " << ksi << std::endl;     
  summary << std::endl;	

  parameters_file << n_ele << " " << ntim << " " << X << " " << a_0 << " " << nu[0] << " "
		  << nu[1] << " " << cs_t << " " << cs_t/ksi << " " << 180*psi/PI << " " 
		  << 180*phi/PI << " ";

  fracture_law->printSelf(parameters_file, summary);
  contact_law->printSelf(parameters_file, summary);}
