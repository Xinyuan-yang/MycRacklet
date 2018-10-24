/* -------------------------------------------------------------------------- */
#include "regularized_coulomb_law.hh"
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it) {


  if (it == 1) {

    sigma_np1[i] = norm_comp_stress;

  }

  else if ((contact_history[i]+1) == it) {

   sigma_np1[i] = (sigma_np1[i] + dt_tstar*norm_comp_stress) / (1 + dt_tstar);

   //norm_comp_stress = sigma_np1[i];

  }

  else {

    sigma_np1[i] = dt_tstar*norm_comp_stress/(1 + dt_tstar);

    //norm_comp_stress = sigma_np1[i];
  }

  strength = cf*fabs(sigma_np1[i]);
  contact_history[i] = it;

}

/* -------------------------------------------------------------------------- */
void RegularizedCoulombLaw::restart(bool pausing, UInt nele_2d) {

  restartData(sigma_np1, "restart_sigma_np1.cra",pausing, nele_2d);
  restartData(contact_history, "restart_contact_history.cra",pausing, nele_2d);
}

