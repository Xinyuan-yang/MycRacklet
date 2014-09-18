/* -------------------------------------------------------------------------- */
#include <vector>
#include <stdio.h>
#include "crack_profile.hh"

/* -------------------------------------------------------------------------- */
template<class Bed>
inline CircularBuffer<Bed>::CircularBuffer(){

}

/* -------------------------------------------------------------------------- */
template<class Bed>
inline CircularBuffer<Bed>::CircularBuffer(std::vector<int> cut){

  resize(cut);

}
/* -------------------------------------------------------------------------- */
template<class Bed>
inline CircularBuffer<Bed>::~CircularBuffer(){

}

/* -------------------------------------------------------------------------- */
template<class Bed>
void CircularBuffer<Bed>::recall(std::vector<std::complex<double> > & res, int it) {
  
  for (int j = 0; j < nele; ++j) {
    
    int cutt = box[2*j].size();
    int i_buff = it - (int)(it/cutt) * cutt;

    res[j] = {box[2*j][i_buff], box[2*j+1][i_buff]};
  }

}

/* -------------------------------------------------------------------------- */
template<class Bed>
void CircularBuffer<Bed>::convolution(std::vector<std::vector<double> > & K, std::vector<std::complex<double> > & res, int it){

  double temp1, temp2;
  int tcut;
  int i_buff;

  for (int j = 0; j < nele; ++j) {
	
    tcut = box[2*j].size();
    temp1 = 0;
    temp2 = 0;

    for (int m = std::max(it - tcut, 0); m < it; ++m) {
	  
      i_buff = m - (int)(m/tcut) * tcut;

      temp1 = K[j][it-m-1] * box[2*j][i_buff]  + temp1;
      temp2 =  K[j][it-m-1] * box[2*j+1][i_buff]  + temp2;
      
    }
	 
    res[j] = {temp1 ,temp2};    
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
void CircularBuffer<Bed>::store(Bed * datas, int it) {

  for (int j = 0; j < nele; ++j) {

    int cutt = box[2*j].size();
    int i_buff = it - (int)(it/cutt) * cutt;

    for (int img = 0; img < 2; ++img) {
  
      box[2*j+img][i_buff] = datas[2*j+img];
    
    } 
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
void CircularBuffer<Bed>::resize(std::vector<int> cut) {

  nele = cut.size();

  box.resize(2*nele);

  for (int i = 0; i < nele; ++i) {
    for (int img= 0; img < 2; ++img) {

      box[2*i+img].resize(cut[i]);
    }
  }
}

/* -------------------------------------------------------------------------- */
template<class Bed>
void CircularBuffer<Bed>::resize(int buf_s, int bed_s) {

  nele = buf_s;  

  box.resize(buf_s);

 for (int i = 0; i < buf_s; ++i) {
   
   box[i].resize(bed_s);

 }

}

