/**
 * @file   spectral_model_inline_impl.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Fri Nov  9 19:10:51 2012
 *
 * @brief  Implementation of the inline functions for the SpectralModel class
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
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(){

  ntim = 0;
  nu.resize(2);
  mu.resize(2);
  rho.resize(2);
  cs.resize(2);
  t_cut.resize(2);
  ksi = 0;
  interface_law=NULL;
}

/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(UInt nele_x, UInt nb_time_steps,
				    Real dom_size_x, Real nu,
				    Real E, Real cs, UInt tcut, 
				    const std::string & simulation_summary,
				    const std::string output_dir) :
  SpectralModel(std::vector<UInt> {nele_x,1},
		nb_time_steps,std::vector<Real> {dom_size_x,0.},
		nu,nu,E,E,cs,cs,tcut,tcut,
		simulation_summary,output_dir){
}

/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(UInt nele_x, UInt nb_time_steps,
				    Real dom_size_x,
				    Real nu_top, Real nu_bot, 
				    Real E_top, Real E_bot,
				    Real cs_top, Real cs_bot,
				    UInt tcut_top, UInt tcut_bot, 
				    const std::string & simulation_summary,
				    const std::string output_dir) :
  SpectralModel(std::vector<UInt> {nele_x,1},
		nb_time_steps,std::vector<Real> {dom_size_x,0.},
		nu_top,nu_bot,E_top,E_bot,
		cs_top,cs_bot,tcut_top,tcut_bot,
		simulation_summary,output_dir){
  
}

/* -------------------------------------------------------------------------- */

inline SpectralModel::SpectralModel(std::vector<UInt> nele, UInt nb_time_steps, 
				    std::vector<Real> dom_size,
				    Real nu, Real E,
				    Real cs, UInt tcut, 
				    const std::string & simulation_summary,
				    const std::string output_dir) :
  SpectralModel(nele,nb_time_steps,
		dom_size,
		nu,nu,E,E,cs,cs,tcut,tcut,
		simulation_summary,output_dir){
}

/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(std::vector<UInt> nele, UInt nb_time_steps, 
				    std::vector<Real> dom_size,
				    Real nu_top, Real nu_bot, 
				    Real E_top, Real E_bot,
				    Real cs_top, Real cs_bot,
				    UInt tcut_top, UInt tcut_bot, 
				    const std::string & simulation_summary,
				    const std::string output_dir) {

  n_ele.resize(2);
  n_ele = nele;
  ntim = nb_time_steps;
  X = dom_size;
  if(n_ele[1]==1){interface_dim=1;}
  else{interface_dim=2;}
  nu.resize(2);
  nu[0] = nu_top;
  nu[1] = nu_bot;
  mu.resize(2);
  mu[0] = E_top/(2*(1+nu_top));
  mu[1] = E_bot/(2*(1+nu_bot));
  rho.resize(2);
  rho[0] = mu[0]/(cs_top*cs_top);
  rho[1] = mu[1]/(cs_bot*cs_bot);
  ksi = cs_top/cs_bot;
  cs.resize(2);
  cs[0] = cs_top;
  cs[1] = cs_bot;
  t_cut.resize(2);
  t_cut[0] = tcut_top;
  t_cut[1] = tcut_bot;
  dx.resize(2);
  q0.resize(2);
  nele_fft.resize(2);
  this->data_initialize(output_dir, simulation_summary);
  interface_law=NULL;
}

/* -------------------------------------------------------------------------- */

template<>
inline void SpectralModel::initFrequency<1>() {

  q00.resize(total_nele_fft-1);
  
  for (UInt k = 0; k < nele_fft[0]-1; ++k) {
    q00[k]=(Real)(k+1)*q0[0];
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void SpectralModel::initFrequency<2>() {

  q00.resize(total_nele_fft-1);
  kk.resize(total_nele_fft-1);
  mm.resize(total_nele_fft-1);

  for (UInt m = 0; m < nele_fft[1]; ++m) {
    for (UInt k = 0; k < nele_fft[0]; ++k) {
      if (k==0 && m==0){continue;}

      if (m<=nele_fft[1]/2) {
      q00[k+nele_fft[0]*m-1] = sqrt((Real)k*q0[0]*(Real)k*q0[0]+(Real)m*q0[1]*(Real)m*q0[1]);
      kk[k+nele_fft[0]*m-1]=(Real)k*q0[0];
      mm[k+nele_fft[0]*m-1]=(Real)m*q0[1];
      }
      else {
	q00[k+nele_fft[0]*m-1] = sqrt((Real)k*q0[0]*(Real)k*q0[0]+(Real)(nele_fft[1]-m)*q0[1]*(Real)(nele_fft[1]-m)*q0[1]);
	kk[k+nele_fft[0]*m-1]=(Real)k*q0[0];
	mm[k+nele_fft[0]*m-1]=-(Real)(nele_fft[1]-m)*q0[1];
      }
    }
  }

  h33u1 = new Real[2*(total_nele_fft-1)];
  h11u3 = new Real[2*(total_nele_fft-1)];
  h12u3 = new Real[2*(total_nele_fft-1)];
}
/* -------------------------------------------------------------------------- */
template<>
inline void SpectralModel::computeConvolutions<1>(Real * F_it, UInt side) {
   
  Real * U_it;
  SpectralConvolutionManager * convo_manager;

  if(side==0) {
    U_it = U_top;
    convo_manager = convo_manager_top;}
  else {
    U_it = U_bot;
    convo_manager = convo_manager_bot;}

  convo_manager->computeConvolution(h11u1,0,0);
  convo_manager->computeConvolution(h12u1,0,2);
  convo_manager->computeConvolution(h12u2,1,2);
  convo_manager->computeConvolution(h22u2,1,1);
  convo_manager->computeConvolution(h33u3,2,3);
 
  Real * t11u1(h11u1);
  Real * t12u1(h12u1);
  Real * t12u2(h12u2);
  Real * t22u2(h22u2);
  Real * t33u3(h33u3);

  Real * U1_it = U_it;
  Real * U2_it = U_it+2*(total_nele_fft-1);
  Real * U3_it = U_it+4*(total_nele_fft-1);

  std::vector<Real> sgn = {1.0, -1.0};

  for (UInt mode = 0; mode < total_nele_fft-1; ++mode) {
    for(UInt i = 0; i < 2; ++i) { //i is related to complex plane

      *(F_it+2+i) = -sgn[side]* (mu[side]*q00[mode])* *(t11u1+i)
	+ sgn[1-i] *(mu[side]*q00[mode])* *(t12u2+1-i)  
	+ sgn[1-i] *(2-eta[side])*mu[side]*(q00[mode])* *(U2_it+1-i);
      
      *(F_it+2*total_nele_fft+2+i) = -sgn[side] * (mu[side]*q00[mode])* *(t22u2+i)
	+sgn[i]* (mu[side]*q00[mode])* *(t12u1+1-i)
	+sgn[i]* (2-eta[side])*mu[side]*(q00[mode]* *(U1_it+1-i));

      *(F_it+4*total_nele_fft+2+i) = 
	- sgn[side]* (mu[side]*q00[mode])* *(t33u3+i); 
    }   
    
    F_it+=2;
    U1_it+=2;
    U2_it+=2;
    U3_it+=2;
    t11u1+=2;	
    t22u2+=2;
    t33u3+=2;
    t12u1+=2;
    t12u2+=2;
  }
}
/* -------------------------------------------------------------------------- */
template<>
inline void SpectralModel::computeConvolutions<2>(Real * F_it, UInt side) {
   
  Real * U_it;
  SpectralConvolutionManager * convo_manager;

  if(side==0) {
    U_it = U_top;
    convo_manager = convo_manager_top;}
  else {
    U_it = U_bot;
    convo_manager = convo_manager_bot;}

  convo_manager->computeConvolution(h11u1,0,0);
  convo_manager->computeConvolution(h12u1,0,2);
  convo_manager->computeConvolution(h12u2,1,2);
  convo_manager->computeConvolution(h22u2,1,1);
  convo_manager->computeConvolution(h33u3,2,3);
  convo_manager->computeConvolution(h11u3,2,0);
  convo_manager->computeConvolution(h12u3,2,2);
  convo_manager->computeConvolution(h33u1,0,3);

  Real * t11u1(h11u1);
  Real * t11u3(h11u3);
  Real * t12u1(h12u1);
  Real * t12u2(h12u2);
  Real * t12u3(h12u3);
  Real * t22u2(h22u2);
  Real * t33u1(h33u1);
  Real * t33u3(h33u3);

  Real * U1_it = U_it;
  Real * U2_it = U_it+2*(total_nele_fft-1);
  Real * U3_it = U_it+4*(total_nele_fft-1);

  std::vector<Real> sgn = {1.0, -1.0};

  for (UInt mode = 0; mode < total_nele_fft-1; ++mode) {
    for(UInt i = 0; i < 2; ++i) { //i is related to complex plane

      *(F_it+2+i) = -sgn[side]* (mu[side]*q00[mode])* *(t11u1+i)*kk[mode]*kk[mode]/(q00[mode]*q00[mode])
	- sgn[side]* (mu[side]*q00[mode])* *(t11u3+i)*kk[mode]*mm[mode]/(q00[mode]*q00[mode]) 
	- sgn[side]* (mu[side]*q00[mode])* *(t33u1+i)*mm[mode]*mm[mode]/(q00[mode]*q00[mode]) 
	+ sgn[side]* (mu[side]*q00[mode])* *(t33u3+i)*(kk[mode]*mm[mode])/(q00[mode]*q00[mode])
	+ sgn[1-i] *(mu[side]*q00[mode])* *(t12u2+1-i)*(kk[mode])/q00[mode]  
	+ sgn[1-i] *(2-eta[side])*mu[side]*(kk[mode])* *(U2_it+1-i);
      
      *(F_it+2*total_nele_fft+2+i) = -sgn[side] * (mu[side]*q00[mode])* *(t22u2+i)
	+sgn[i]* (mu[side]*q00[mode])* *(t12u1+1-i)*kk[mode]/q00[mode]
	+sgn[i]* (mu[side]*q00[mode])* *(t12u3+1-i)*mm[mode]/q00[mode]
	+sgn[i]* (2-eta[side])*mu[side]*(kk[mode]* *(U1_it+1-i)+mm[mode]* *(U3_it+1-i));

      *(F_it+4*total_nele_fft+2+i) = 
	- sgn[side]* (mu[side]*q00[mode])* *(t11u1+i)*kk[mode]*mm[mode]/(q00[mode]*q00[mode]) 
	- sgn[side]* (mu[side]*q00[mode])* *(t11u3+i)*mm[mode]*mm[mode]/(q00[mode]*q00[mode]) 
	+ sgn[side]* (mu[side]*q00[mode])* *(t33u1+i)*(kk[mode]*mm[mode])/(q00[mode]*q00[mode]) 
	- sgn[side]* (mu[side]*q00[mode])* *(t33u3+i)*kk[mode]*kk[mode]/(q00[mode]*q00[mode]) 
	+sgn[1-i]*(mu[side]*q00[mode])* *(t12u2+1-i)*(mm[mode])/q00[mode]  
	+sgn[1-i]*(2-eta[side])*mu[side]*(mm[mode])* *(U2_it+1-i);
    }   

    F_it+=2;
    U1_it+=2;
    U2_it+=2;
    U3_it+=2;
    t11u1+=2;	
    t22u2+=2;
    t33u3+=2;
    t12u1+=2;
    t12u2+=2;
    t33u1+=2;
    t11u3+=2;
    t12u3+=2;
  }
}
/* -------------------------------------------------------------------------- */
inline SpectralModel::~SpectralModel(){

  for (UInt side = 0; side < 2; ++side) {
    displacements[side].finalizeFFT();
    stresses[side].finalizeFFT();
  }

  this->data_finalize();
  
  delete displ_jump;
  delete veloc_jump;

  if(h11u1){
    delete[] h11u1;
    delete[] h22u2;
    delete[] h33u3;
    delete[] h12u1;
    delete[] h12u2;
    if(interface_dim==2) {
      delete[] h12u3;
      delete[] h33u1;
      delete[] h11u3;
    }
    delete[] U_top;
    delete[] U_bot;
    delete[] F_k;
  }

  delete convo_manager_top;
  delete convo_manager_bot;

};
