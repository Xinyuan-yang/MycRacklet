/* -------------------------------------------------------------------------- */
#include "coulomb_law.hh"
#include <stdlib.h>
#include <math.h>
#include <fstream>
/* -------------------------------------------------------------------------- */
void CoulombLaw::computeFricStrength(double & norm_comp_stress, double & strength, int i, int it) {

  strength = cf*fabs(norm_comp_stress);

}

/* -------------------------------------------------------------------------- */
void CoulombLaw::printSelf(std::ofstream & parameters_file, std::ofstream & summary) {

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " CONTACT LAW VARIABLES " << std::endl;
  summary << "* Type of contact law: Coulomb law" << std::endl;		        
  summary << "* Coefficient of friction: " << cf << std::endl;
  summary << std::endl;	

  parameters_file << cf << " ";
}
