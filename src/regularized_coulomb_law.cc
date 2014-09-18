/* -------------------------------------------------------------------------- */
#include "regularized_coulomb_law.hh"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::computeFricStrength(double & norm_comp_stress, double & strength, int i, int it) {


  if (it == 1) {

    sigma_np1[i] = norm_comp_stress;

  }

  else if ((contact_history[i]+1) == it) {

    sigma_np1[i] = (sigma_np1[i] + dt_tstar*norm_comp_stress) / (1 + dt_tstar);

    norm_comp_stress = sigma_np1[i];

  }

  else {
    
    sigma_np1[i] = dt_tstar*norm_comp_stress/(1 + dt_tstar);

    norm_comp_stress = sigma_np1[i];
  }

  strength = cf*fabs(norm_comp_stress);
  contact_history[i] = it;

}

/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::printSelf(std::ofstream & parameters_file, std::ofstream & summary) {

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " CONTACT LAW VARIABLES " << std::endl;
  summary << "* Type of contact law: Coulomb law with normal pressure regularization" << std::endl;	        
  summary << "* Coefficient of friction: " << cf << std::endl;
  summary << "* Regularization time scale dt/t*:  " << dt_tstar << std::endl; 
  summary << std::endl;	

  parameters_file << cf << " ";

}
