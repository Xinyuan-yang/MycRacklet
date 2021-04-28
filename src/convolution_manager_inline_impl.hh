/* -------------------------------------------------------------------------- */
template<class Bed>
void ConvolutionManager::KernelFunctor::readBinary(std::ifstream & file, Bed & val){

char * val_char = reinterpret_cast<char *>(&val);
  file.read((char*)val_char, sizeof(Bed));
}

/* -------------------------------------------------------------------------- */
template<>
inline void ConvolutionManager::KernelFunctor::readBinary<std::vector<Real> >(std::ifstream & file, 
										std::vector<Real> & val){
  
  UInt numt = val.size();
  Real * tmp_values = new Real[2*numt];
  
  file.read((char*) tmp_values, 2*numt*sizeof(Real));
  for (UInt i = 0; i < numt; ++i) {
    tmp_values[i] = tmp_values[2*i+1];
  }
  
  for (UInt i = 0; i < numt; ++i) {
    val[i] = tmp_values[i];
  }

  delete[] tmp_values;
}

 
/* -------------------------------------------------------------------------- */
inline void ConvolutionManager::preintegrateKernel(UInt i, KernelFunctor * funct, 
						   UInt k_start, Real alpha) {

  Real tau = (Real)i*dx*alpha;
  Real k_tau = (*funct)(tau);
  Real next = ((Real)i+1)*dx*alpha;
  Real k_next = (*funct)(next);

  *(K+k_start+i) = 0.5*(k_tau + k_next)*dx*alpha;
}

/* -------------------------------------------------------------------------- */
inline void ConvolutionManager::storeFields(Real new_val) {

  UInt n = field->getStep();
  if ((n<size)&&(func))
    preintegrateKernel(n,this->func);  
  storeFields(*field, new_val);

}

/* -------------------------------------------------------------------------- */
inline void ConvolutionManager::storeFields(RingBuffer<Real> & buffer_destination,
					    Real new_val) {
  buffer_destination << new_val;
}

/* -------------------------------------------------------------------------- */
inline Real ConvolutionManager::computeConvolution() {

  return convolute(*field,size);
}

/* -------------------------------------------------------------------------- */
inline Real ConvolutionManager::convolute(RingBuffer<Real> & buffer, 
					    UInt cut, UInt k_start) {

  UInt n = buffer.getStep();
  UInt max_n = std::min(n,cut);
  Real convo = 0.;

  Real * it;
  Real * it_K = K+k_start+max_n-1;
  Real * current;
  Real * begin;
  Real * end;

  current = buffer.current()+1;
  begin = buffer.begin();
  end = buffer.end()+1;

  for(it = current; it!=end; ++it) {

    convo += *it * *it_K;
    --it_K;
  }
  
  for(it = begin; it!=current; ++it){

    convo += *it * *it_K;
    --it_K;
  }
  return convo;
}
