/**
 * @file   convolution_manager.cc
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Fri Nov 13 12:15:34 2015
 *
 * @brief  Implementation of the ConvolutionManager class
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
#include "convolution_manager.hh"
#include <cmath>
#include <assert.h>
#include <limits>
/* -------------------------------------------------------------------------- */
Real ConvolutionManager::KernelFunctor::operator()(Real tau) {

  UInt i_low;
  Real val_low;
  Real val_up;
  Real ker_low;
  Real ker_up;
  Real result;

  i_low = (UInt)(tau/delta_smpl);

  ker_low = kernel[i_low];
  ker_up = kernel[i_low+1];

  val_low = ((Real) i_low) * delta_smpl;
  val_up = ((Real) i_low +1)*delta_smpl;

  result = (ker_up-ker_low) * (tau-val_low)/(val_up-val_low) + ker_low;

  return result;

}
/* -------------------------------------------------------------------------- */
Real ConvolutionManager::KernelFunctor::loadKernelFromFile(std::string filename) {
  
  std::ifstream file (filename, std::ios::binary);

  if (!file.is_open()){
    std::stringstream err;
    err << "Unable to open kernel file" << std::endl;
    err << "Check that the file " << filename
	<< " is in the current folder" << std::endl;
    cRacklet::error(err);
  }
  
  Real cdcs_f;
  Real nu_f;
  UInt numt;
  Real kernel_cut;

  readBinary(file, numt);
  readBinary(file, kernel_cut);

  readBinary(file, numt);
  readBinary(file, delta_smpl);
  readBinary(file, cdcs_f);
  readBinary(file, nu_f);
 
  kernel.resize(numt);
  readBinary(file, kernel);
  
  for (UInt i = 0; i < numt; ++i) {

    if (std::isnan(kernel[i])) {
      std::stringstream err;
      err << "!! NaN values read for kernel at the " << i << "th values !!" << std::endl;
      cRacklet::error(err);
    }  
  }
  
  return nu_f;
}

/* -------------------------------------------------------------------------- */
ConvolutionManager::~ConvolutionManager() {

  if(this->func) {
    delete this->func;
    delete this->field;
  }
  if(this->K)
    delete[] K;
  if(this->field_values)
    delete[] field_values;
}

/* -------------------------------------------------------------------------- */
void ConvolutionManager::init(UInt cut) {

  allocateMemory(cut);
  this->field = new RingBuffer<Real>;
  initRingBuffer(*field, field_values, cut);
}
/* -------------------------------------------------------------------------- */
void ConvolutionManager::allocateMemory(Idx size) {

  this->size = size;
  field_values = new Real[size];
}

/* -------------------------------------------------------------------------- */
void ConvolutionManager::initRingBuffer(RingBuffer<Real> & buffer_to_init,
					Real * buffer_start, UInt buffer_length) {

  buffer_to_init.init(buffer_start, buffer_length, 0.);
}

/* -------------------------------------------------------------------------- */
void ConvolutionManager::initK() {

    K = new Real[size];
}
/* -------------------------------------------------------------------------- */
void ConvolutionManager::loadKernel(KernelFunctor * func) {
  this->func = func;
  initK();
}

/* -------------------------------------------------------------------------- */
Real ConvolutionManager::loadKernel(std::string filename) {
  
  this->func=new KernelFunctor();
  initK();
  return func->loadKernelFromFile(filename);
}

/* -------------------------------------------------------------------------- */
void ConvolutionManager::dumpKernel(std::string filename, Real from, Real to, UInt nb_points) {

  std::ofstream outfile; 
  outfile.open(filename.c_str());
  Real dtau = (to-from)/(Real)(nb_points);
  Real tau = from;

  while (tau<to) {

    outfile << tau << " " << (*func)(tau) << std::endl;
    tau += dtau;    
  }
}
