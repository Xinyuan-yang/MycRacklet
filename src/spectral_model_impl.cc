/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(){

  n_ele = 0;
  X = 0;
  ntim = 0;
  nu.resize(2);
  mu.resize(2);
  t_cut.resize(2);
  ksi = 0;
  overlapping = 0;
}

/* -------------------------------------------------------------------------- */
inline SpectralModel::SpectralModel(int nele, int nb_time_steps, double dom_size, 
				    double crack_size, double nu_top, double nu_bot, 
				    double E_top, double E_bot, double cs_top, 
				    double cs_bot, int tcut_top, int tcut_bot, 
				    bool overlap, unsigned int l_index, 
				    FractureLaw * fracturelaw, ContactLaw * contactlaw) {


  n_ele = nele;
  ntim = nb_time_steps;
  X = dom_size;
  a_0 = crack_size;
  nu.resize(2);
  nu[0] = nu_top;
  nu[1] = nu_bot;
  mu.resize(2);
  mu[0] = E_top/(2*(1+nu_top));
  mu[1] = E_bot/(2*(1+nu_bot));
  ksi = cs_top/cs_bot;
  cs_t = cs_top;
  t_cut.resize(2);
  t_cut[0] = tcut_top;
  t_cut[1] = tcut_bot;
  overlapping = overlap;
  load_index = l_index;
  fracture_law = fracturelaw;
  contact_law = contactlaw;

  
} 

/* -------------------------------------------------------------------------- */
inline SpectralModel::~SpectralModel(){

  delete displ_jump;
  delete veloc_jump;
};


/* -------------------------------------------------------------------------- */
template<class Bed>
void SpectralModel::readBinary(std::ifstream & file, Bed & val){

char * val_char = reinterpret_cast<char *>(&val);
  file.read((char*)val_char, sizeof(Bed));
}

/* -------------------------------------------------------------------------- */
template<>
inline void SpectralModel::readBinary<std::vector<double> >(std::ifstream & file, 
				      std::vector<double> & val){
  
  int numt = val.size();
  double * tmp_values = new double[2*numt];

  file.read((char*) tmp_values, 2*numt*sizeof(double));
  for (int i = 0; i < numt; ++i) {
    tmp_values[i] = tmp_values[2*i+1];
  }
  
  for (int i = 0; i < numt; ++i) {
    val[i] = tmp_values[i];
  }

  delete[] tmp_values;
}
