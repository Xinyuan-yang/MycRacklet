/* -------------------------------------------------------------------------- */
#include "cohesive_law.hh"
#include <fstream>
/* -------------------------------------------------------------------------- */
void CohesiveLaw::updateFractureLaw(std::vector<double> & nor_strength, std::vector<double> & shr_strength,
				    std::vector<unsigned int> & ind_crack, CrackProfile & nor_opening, 
				    CrackProfile & shr_opening) {


  int n_ele = nor_strength.size();
  double aux;
  double max_n_str;
  double max_s_str;

  for (int i = 0; i < n_ele; ++i) {
  

    if (ind_crack[i] == 4) {

      nor_strength[i] = 1000*max_nor_strength;
      shr_strength[i] = 1000*max_shr_strength;

    }
    
    else {

      	max_n_str = max_nor_strength;
	max_s_str = max_shr_strength;
      
      }

    if ((nor_opening[i]==0)&&(shr_opening[i]==0)){

      nor_strength[i] = max_n_str;
      shr_strength[i] = max_s_str;

    }
    
    else {
      
      aux = sqrt((nor_opening[i]/crit_nor_opening)*(nor_opening[i]/crit_nor_opening)+
		 (shr_opening[i]/crit_shr_opening)*(shr_opening[i]/crit_shr_opening));
      
      if ((aux>=1)||(nor_strength[i] * shr_strength[i] == 0)) {

	if ((nor_strength[i] == 0)&&(nor_opening[i] <=0)) ind_crack[i] = 3;
	
	else ind_crack[i] = 2;
	
	nor_strength[i] = 0.0;
	shr_strength[i] = 0.0;
	
	}

      else {

	  nor_strength[i] = max_n_str * (1-aux);
	  shr_strength[i] = max_s_str * (1-aux);
	  ind_crack[i] = 1;
      }
    }
  }    
}

/* -------------------------------------------------------------------------- */
void CohesiveLaw::printSelf(std::ofstream & parameters_file, std::ofstream & summary) {

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " FRACTURE LAW VARIABLES " << std::endl;
  summary << "* Type of fracture law: Linear rate-independant coupled cohesive law" 
	  << std::endl;		        
  summary << "* Maximum shear strength: " << max_shr_strength << std::endl;
  summary << "* Maximum normal strength: " << max_nor_strength << std::endl;
  summary << "* Critical shear opening: " << crit_shr_opening << std::endl;
  summary << "* Critical normal opening: " << crit_nor_opening << std::endl;
  summary << std::endl;	

  parameters_file << max_nor_strength << " " << crit_nor_opening << " " << max_shr_strength 
		  << " " << crit_shr_opening << " ";
  
}

