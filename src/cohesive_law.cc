/* -------------------------------------------------------------------------- */
#include "cohesive_law.hh"
#include <fstream>
#include <algorithm>
/* -------------------------------------------------------------------------- */
void CohesiveLaw::updateFractureLaw(std::vector<Real> & nor_strength, std::vector<Real> & shr_strength,
				    std::vector<UInt> & ind_crack, CrackProfile & nor_opening, 
				    CrackProfile & shr_opening) {


  UInt n_ele = max_nor_strength.size();
  Real aux;

  for (UInt i = 0; i < n_ele; ++i) {     

    if ((nor_opening[i]==0)&&(shr_opening[i]==0)){

      nor_strength[i] = max_nor_strength[i];
      shr_strength[i] = max_shr_strength[i];

    }
    
    else {
      
      aux = sqrt((nor_opening[i]/crit_nor_opening[i])*(nor_opening[i]/crit_nor_opening[i])+
		 (shr_opening[i]/crit_shr_opening[i])*(shr_opening[i]/crit_shr_opening[i]));
      
      if ((aux>=1)||(nor_strength[i] * shr_strength[i] == 0)) {

	if ((nor_strength[i] == 0)&&(cRacklet::is_negative(nor_opening[i]))){} 
	//nodes whith contact are handled by ContactLaw} 
	
	else 
	  ind_crack[i] = 2;
	
	nor_strength[i] = 0.0;
	shr_strength[i] = 0.0;
	
      }

      else {
	
	nor_strength[i] = max_nor_strength[i] * (1-aux);
	shr_strength[i] = max_shr_strength[i] * (1-aux);
	ind_crack[i] = 1;

      }
    }
  }    
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::printSelf(std::ofstream & parameters_file, std::ofstream & summary) {

  Real max_nor_str = *std::max_element(max_nor_strength.begin(),max_nor_strength.end());
  Real min_nor_str = *std::min_element(max_nor_strength.begin(),max_nor_strength.end());

  Real max_shr_str = *std::max_element(max_shr_strength.begin(),max_shr_strength.end());
  Real min_shr_str = *std::min_element(max_shr_strength.begin(),max_shr_strength.end());

  Real max_nor_op = *std::max_element(crit_nor_opening.begin(),crit_nor_opening.end());
  Real min_nor_op = *std::min_element(crit_nor_opening.begin(),crit_nor_opening.end());

  Real max_shr_op = *std::max_element(crit_shr_opening.begin(),crit_shr_opening.end());
  Real min_shr_op = *std::min_element(crit_shr_opening.begin(),crit_shr_opening.end());

  bool homog = false;

  if ((max_nor_str==min_nor_str)&&(max_shr_str==min_shr_str)
      &&(max_nor_op==min_nor_op)&&(max_shr_str==min_shr_str)) homog = true;

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " FRACTURE LAW VARIABLES " << std::endl;
  summary << "* Type of fracture law: Linear rate-independant coupled cohesive law" 
	  << std::endl;		        
  summary << "* Maximum shear strength: " << min_shr_str << std::endl;
  summary << "* Maximum normal strength: " << min_nor_str << std::endl;
  summary << "* Critical shear opening: " << min_shr_op << std::endl;
  summary << "* Critical normal opening: " << min_nor_op << std::endl;
  summary << std::endl;	
  parameters_file << min_nor_str << " " << min_nor_op << " " << min_shr_str 
		  << " " << min_shr_op << " ";
  if (!homog) {
  summary << "* Heterogeneous tougher area with properties: " << std::endl;  
  summary << "* Maximum shear strength: " << max_shr_str << std::endl;
  summary << "* Maximum normal strength: " << max_nor_str << std::endl;
  summary << "* Critical shear opening: " << max_shr_op << std::endl;
  summary << "* Critical normal opening: " << max_nor_op << std::endl;
  summary << std::endl;
  parameters_file << max_nor_str << " " << max_nor_op << " " << max_shr_str 
		  << " " << max_shr_op << " ";
  }
}

