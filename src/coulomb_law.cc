/* -------------------------------------------------------------------------- */
#include "coulomb_law.hh"
#include <stdlib.h>
#include <math.h>
#include <fstream>
/* -------------------------------------------------------------------------- */
void CoulombLaw::computeFricStrength(Real & norm_comp_stress, Real & strength, UInt i, UInt it) {

  strength = cf*fabs(norm_comp_stress);

}

