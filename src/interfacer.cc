/* -------------------------------------------------------------------------- */
#include "interfacer.hh"
/* -------------------------------------------------------------------------- */
void Interfacer::createCenteredCrack(double max_nor_strength, double max_shr_strength){

  double pos = 0.5*dx;
  double crk_srt = 0.5*(n_ele*dx-crack_size);
  
  for (int i = 0; i < n_ele; ++i) {

    if ((pos <= (crk_srt+crack_size)) && (pos >= (crk_srt))){
      nor_strength[i]=0;
      shr_strength[i]=0;
      ind_crack[i] = 2;
    }
    else {
      nor_strength[i]=max_nor_strength;
      shr_strength[i]=max_shr_strength;
      ind_crack[i] = 0;
    }
    pos +=dx;
  }
}


/* -------------------------------------------------------------------------- */
void Interfacer::createLeftSidedCrack(double max_nor_strength, double max_shr_strength){

  double pos = 1.5*dx;

  for (int i = 0; i < n_ele; ++i) {
    
    if (pos <= crack_size){
      nor_strength[i]=0;
      shr_strength[i]=0;
      ind_crack[i] = 2;
    }

    else if (pos > 0.9*n_ele*dx) {

      nor_strength[i]=1000*max_nor_strength;
      shr_strength[i]=1000*max_shr_strength;
      ind_crack[i]=4;
    }

    else {
      nor_strength[i]=max_nor_strength;
      shr_strength[i]=max_shr_strength;
      ind_crack[i] = 0;
    }
    pos +=dx;
  }
}

/* -------------------------------------------------------------------------- */
void Interfacer::createIncohIntfc() {

  for (int i = 0; i < n_ele; ++i) {
    
      nor_strength[i]=0;
      shr_strength[i]=0;
      ind_crack[i] = 2;
  }
}
