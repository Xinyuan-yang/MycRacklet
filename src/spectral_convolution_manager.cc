/* -------------------------------------------------------------------------- */
#include "spectral_convolution_manager.hh"
#if defined (_OPENMP)
#include <omp.h>
#endif
/* -------------------------------------------------------------------------- */
SpectralConvolutionManager::SpectralConvolutionManager(Real dx, std::vector<Real> & modal_number,
						       UInt nb_fields, UInt nb_kernels)
  : ConvolutionManager(dx) {
  this->nb_modes = modal_number.size();
  this->nb_fields = nb_fields;
  this->nb_kernels = nb_kernels;
  this->mod_cut.resize(nb_modes);
  this->field = new RingBuffer<Real>[2*nb_modes*nb_fields];
  this->functors.resize(nb_kernels);
  this->k_start.resize(nb_modes);
  this->j_ksi.resize(nb_modes);
  for (UInt i=0; i < nb_modes; ++i) {
    j_ksi[i] = modal_number[i];
  }
}

/* -------------------------------------------------------------------------- */
SpectralConvolutionManager::~SpectralConvolutionManager() {

  delete[] field;
  for (UInt k = 0; k < nb_kernels; ++k) {
    if(functors[k])
      delete functors[k];
  }
}
/* -------------------------------------------------------------------------- */
void SpectralConvolutionManager::init(Real cut, UInt nb_time) {

  this->size_per_field = 0;
  
  if(nb_time==0) //In case of unknown number of time steps
    nb_time=std::numeric_limits<UInt>::max();

  for (UInt i=0; i < nb_modes; ++i) {
    mod_cut[i] = std::min(nb_time,(UInt)(cut/(dx*j_ksi[i]))-1);
    k_start[i] = size_per_field;
    size_per_field += mod_cut[i];
  }

  Idx total_size = 2*size_per_field*nb_fields;

  ConvolutionManager::allocateMemory(2*size_per_field*nb_fields);

  Real * it = field_values;

  for(UInt f=0; f < nb_fields; ++f) {
    for (UInt i=0; i < nb_modes; ++i) {
      for (UInt j=0; j<2; ++j){
	ConvolutionManager::initRingBuffer(*(field+2*i+j+2*nb_modes*f), it, mod_cut[i]);
	it += mod_cut[i];
      }
    }
  }

  initK();
}

/* -------------------------------------------------------------------------- */
ConvolutionManager::KernelFunctor ** SpectralConvolutionManager::getKernelFunctor(UInt k) {

  return &functors[k];
}

/* -------------------------------------------------------------------------- */
Real SpectralConvolutionManager::loadKernel(std::string filename, UInt kernel_id) {

  functors[kernel_id] = new KernelFunctor();
  return functors[kernel_id]->loadKernelFromFile(filename);
}

/* -------------------------------------------------------------------------- */
void SpectralConvolutionManager::initK() {
  
  K = new Real[size_per_field*nb_kernels];
}

/* -------------------------------------------------------------------------- */
void SpectralConvolutionManager::storeFields(Real * new_vals) {

  UInt n = field->getStep();
  Real * it(new_vals);

  for (UInt i=0; i < nb_modes; ++i) {
    if (n<mod_cut[i]) {
      for (UInt k=0; k < nb_kernels; ++k) {
	preintegrateKernel(n,functors[k], k_start[i]+k*size_per_field,j_ksi[i]);
      }
    }
  }
  
  for(UInt f=0; f < nb_fields; ++f) {
    for (UInt i=0; i < nb_modes; ++i) {
      for (UInt j=0; j<2; ++j){
	ConvolutionManager::storeFields(*(field+2*i+j+2*nb_modes*f), *it);
	++it;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void SpectralConvolutionManager::computeConvolution(Real * res, UInt field_id,
						    UInt kernel_id) {

  UInt i,j,chunk,nb_threads;
#if defined (_OPENMP) 
#pragma omp parallel default(shared) private(i,j,chunk,nb_threads)
  {
    nb_threads=omp_get_num_threads();
    chunk=std::max((UInt)1,nb_modes/(50*nb_threads));

#pragma omp for schedule(dynamic,chunk) collapse(2)
#endif
    for (i=0; i < nb_modes; ++i) {
      for (j=0; j<2; ++j){
	*(res+2*i+j) = ConvolutionManager::convolute(*(field+2*i+j+2*nb_modes*field_id), 
					    mod_cut[i], k_start[i]+kernel_id*size_per_field);
      }
    }
#if defined (_OPENMP)
  }
#endif
}

