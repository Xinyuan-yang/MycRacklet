#include "cohesive_law.hh"
#include <random>
#include <chrono>
/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createUniformInterface(Real crit_nor_opening, 
									 Real crit_shr_opening, 
									 Real max_nor_strength, 
									 Real max_shr_strength) {

  std::fill(nor_strength->begin(), nor_strength->end(), max_nor_strength);
  std::fill(shr_strength->begin(), shr_strength->end(), max_shr_strength);
  crit_n_open.resize(total_n_ele, crit_nor_opening);
  crit_s_open.resize(total_n_ele, crit_shr_opening);

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " UNIFORM INTERFACE " << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl;  

  
}

/* -------------------------------------------------------------------------- */
#ifdef CRACKLET_USE_LIBSURFER
template<>
inline void Interfacer<_linear_coupled_cohesive>::createBrownianHeterogInterface(Real crit_nor_opening, 
										 Real crit_shr_opening, 
										 Real max_nor_strength, 
										 Real max_shr_strength,
										 Real rms, long int seed,
										 Real hurst, UInt q0,
										 UInt q1, UInt q2){

  SurfaceGeneratorFilterFFT surf_gen;
  int & grid_size = surf_gen.getGridSize();
  grid_size = std::max(n_ele[0],n_ele[1]);
  Real & Hurst = surf_gen.getHurst();
  Hurst = hurst;
  Real & RMS = surf_gen.getRMS();
  RMS = rms;
  int & Q0 = surf_gen.getQ0();
  Q0 = q0;
  int & Q1 = surf_gen.getQ1();
  Q1 = q1;
  int & Q2 = surf_gen.getQ2();
  Q2 = q2;
  long int & Seed = surf_gen.getRandomSeed();
  Seed = seed;
  surf_gen.Init();
  Surface<Real> & surface = surf_gen.buildSurface();
  std::cout << "Successfully generated surface with an RMS of " << SurfaceStatistics::computeStdev(surface) << std::endl;

  for (UInt ix = 0; ix < n_ele[0]; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {
      
      UInt index = ix+iz*n_ele[0];
      (*nor_strength)[index] = max_nor_strength+surface(ix,iz);      
      (*shr_strength)[index] = max_shr_strength+surface(ix,iz);      
      crit_n_open[index]  = crit_nor_opening;
      crit_s_open[index]  = crit_shr_opening;
    }
  }
}

#endif
/* -------------------------------------------------------------------------- */
template<>inline
void Interfacer<_linear_coupled_cohesive>::createNormalDistributedInterface(Real crit_nor_opening, 
									    Real crit_shr_opening, 
									    Real max_nor_strength, 
									    Real max_shr_strength,
									    Real stddev, Real seed) {
 
  // construct a trivial random generator engine from a time-based seed:
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<Real> distribution (max_shr_strength,stddev);

  for (UInt i = 0; i < total_n_ele; ++i) {

    Real random_strength = distribution(generator);

    (*nor_strength)[i] = random_strength;
    (*shr_strength)[i] = random_strength;
    crit_n_open[i] = crit_nor_opening;
    crit_s_open[i] = crit_shr_opening;
    (*ind_crack)[i] = 0;
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughArea(Real area_start, Real area_end,
								    UInt cracking_index,
								    Real ratio_max_nor_strength, 
								    Real ratio_max_shr_strength,
								    Real ratio_crit_nor_opening, 
								    Real ratio_crit_shr_opening,
								    bool vrtr) {

  UInt i_start = (UInt)(area_start/dx[0]);
  UInt i_end = (UInt)(area_end/dx[0]);

  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt index = ix+iz*n_ele[0];

      (*nor_strength)[index] *= (vrtr+ratio_max_nor_strength);      
      (*shr_strength)[index] *= (vrtr+ratio_max_shr_strength);
      crit_n_open[index] *= (vrtr+ratio_crit_nor_opening);
      crit_s_open[index] *= (vrtr+ratio_crit_shr_opening);
      (*ind_crack)[index] = cracking_index;
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughCrack(Real crack_start, Real crack_end) {

  createThroughArea(crack_start, crack_end, 2, 0., 0.);
  Real initial_crack_size = crack_end-crack_start;
  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Crack starts at: " << crack_start << std::endl
	      << "* Crack ends at: " << crack_end << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "a_l " << crack_start << std::endl
		 << "a_r " << crack_end << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughWall(Real wall_start, Real wall_end) {
  
  createThroughArea(wall_start, wall_end, 6, 1000., 1000.);

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " THROUGH WALL " << std::endl
	      << "* Wall starts at: " << wall_start << std::endl
	      << "* Wall ends at: " << wall_end << std::endl
	      << std::endl;	
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughPolarAsperity(Real position, Real width,
									     Real delta_max_nor_strength, 
									     Real delta_max_shr_strength,
									     Real delta_crit_nor_opening, 
									     Real delta_crit_shr_opening, 
									     bool polarity) {

  if(polarity) {
    createThroughArea(position-width/2,position,5,-delta_max_nor_strength,-delta_max_shr_strength,
		      -delta_crit_nor_opening, -delta_crit_shr_opening, 1);
    createThroughArea(position,position+width/2,4,delta_max_nor_strength,delta_max_shr_strength,
		      delta_crit_nor_opening, delta_crit_shr_opening, 1);
  }
  else {
    createThroughArea(position-width/2,position,4,delta_max_nor_strength,delta_max_shr_strength,
		      delta_crit_nor_opening, delta_crit_shr_opening, 1);
    createThroughArea(position,position+width/2,5,-delta_max_nor_strength,-delta_max_shr_strength,
		      -delta_crit_nor_opening, -delta_crit_shr_opening, 1);
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughMultiPolAsperity(Real start, Real end,
										Real number,
										Real delta_max_nor_strength, 
										Real delta_max_shr_strength,
										Real delta_crit_nor_opening, 
										Real delta_crit_shr_opening, 
										bool polarity) {

  UInt i_start = (UInt)(start/dx[0]);
  UInt i_end = (UInt)(end/dx[0]);

  UInt nb_dx = (UInt)((i_end-i_start)/(2*number));
  if(nb_dx == 0) {
    nb_dx = 1;
    std::cout << "! Asperity width was to small and is currently set to grid size !" << std::endl;
  }
  std::cout << "Asperity width (strong+weak) = " << 2*nb_dx*dx[0] << std::endl;
  Real tough_ratio = (1.+delta_max_shr_strength)/(1.-delta_max_shr_strength);
  std::cout << "Asperity toughness ratio = " << tough_ratio << std::endl;
  
  for (UInt ix = i_start; ix < i_end; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt no = ix-i_start;
      UInt asp_bool = ((no-no%nb_dx)%(2*nb_dx))/nb_dx;
      UInt index = ix+iz*n_ele[0];
      
      if ((asp_bool==1)) {

	(*nor_strength)[index] *= (1+delta_max_nor_strength);      
	(*shr_strength)[index] *= (1+delta_max_shr_strength);
	crit_n_open[index] *= (1+delta_crit_nor_opening);
	crit_s_open[index] *= (1+delta_crit_shr_opening);
	(*ind_crack)[index] = 5;
      }
      else {
	(*nor_strength)[index] *= (1-delta_max_nor_strength);      
	(*shr_strength)[index] *= (1-delta_max_shr_strength);
	crit_n_open[index] *= (1-delta_crit_nor_opening);
	crit_s_open[index] *= (1-delta_crit_shr_opening);
	(*ind_crack)[index] = 4;
      }
    }
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " HETEROGENEOUS INTERFACE " << std::endl
	      << "* Delta critical normal opening: " << delta_crit_nor_opening << std::endl
	      << "* Delta critical shear opening: " << delta_crit_shr_opening << std::endl
	      << "* Delta maximal normal strength: " << delta_max_nor_strength << std::endl
	      << "* Delta maximal shear strength: " << delta_max_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "delta_delta_c_nor " << delta_crit_nor_opening << std::endl
		 << "delta_delta_c_shr " << delta_crit_shr_opening << std::endl
		 << "delta_tau_max_nor " << delta_max_nor_strength << std::endl
		 << "delta_tau_max_shr " << delta_max_shr_strength << std::endl;
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::applyInterfaceCreation() {

  fracture_law = new CohesiveLaw(crit_n_open, crit_s_open, *nor_strength, *shr_strength);
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughCenteredCrack(Real initial_crack_size,
									     Real crit_nor_opening, 
									     Real crit_shr_opening, 
									     Real max_nor_strength, 
									     Real max_shr_strength){
  
  Real pos = 0.5*dx[0];
  Real crk_srt = 0.5*(n_ele[0]*dx[0]-initial_crack_size);
 
  crit_n_open.resize(total_n_ele, crit_nor_opening);
  crit_s_open.resize(total_n_ele, crit_shr_opening);
  
  for (UInt ix = 0; ix < n_ele[0]/2; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt right_cmpnt = ix+n_ele[0]/2+iz*n_ele[0];
      UInt left_cmpnt = n_ele[0]/2-ix-1+iz*n_ele[0];

      if (pos <= initial_crack_size/2){

	(*nor_strength)[left_cmpnt]=0;
	(*shr_strength)[left_cmpnt]=0;
	(*ind_crack)[left_cmpnt] = 2;

	(*nor_strength)[right_cmpnt]=0;
	(*shr_strength)[right_cmpnt]=0;
	(*ind_crack)[right_cmpnt] = 2;
      }
      else {
	(*nor_strength)[left_cmpnt]=max_nor_strength;
	(*shr_strength)[left_cmpnt]=max_shr_strength;
	(*ind_crack)[left_cmpnt] = 0;

	(*nor_strength)[right_cmpnt]=max_nor_strength;
	(*shr_strength)[right_cmpnt]=max_shr_strength;
	(*ind_crack)[right_cmpnt] = 0;
      }
    }
    pos +=dx[0];
  }

  out_summary << "/* -------------------------------------------------------------------------- */ "
	      << std::endl
	      << " CENTERED CRACK " << std::endl
	      << "* Initial crack size: " << initial_crack_size << std::endl
	      << "* Critical normal opening: " << crit_nor_opening << std::endl
	      << "* Critical shear opening: " << crit_shr_opening << std::endl
	      << "* Maximal normal strength: " << max_nor_strength << std::endl
	      << "* Maximal shear strength: " << max_shr_strength << std::endl
	      << std::endl;	

  out_parameters << "a0 " << initial_crack_size << std::endl
		 << "delta_c_nor " << crit_nor_opening << std::endl
		 << "delta_c_shr " << crit_shr_opening << std::endl
		 << "tau_max_nor " << max_nor_strength << std::endl
		 << "tau_max_shr " << max_shr_strength << std::endl;  
}


/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createThroughLeftSidedCrack(Real initial_crack_size,
									      Real crit_nor_opening, 
									      Real crit_shr_opening, 
									      Real max_nor_strength, 
									      Real max_shr_strength){

  Real pos = 0.5*dx[0];

  crit_n_open.resize(total_n_ele, crit_nor_opening);
  crit_s_open.resize(total_n_ele, crit_shr_opening);

  for (UInt ix = 0; ix < n_ele[0]; ++ix) {
    for (UInt iz = 0; iz < n_ele[1]; ++iz) {

      UInt cmpnt = ix+iz*n_ele[0];
    
      if (pos <= initial_crack_size){
	(*nor_strength)[cmpnt]=0;
	(*shr_strength)[cmpnt]=0;
	(*ind_crack)[cmpnt] = 2;
      }
      else {
	(*nor_strength)[cmpnt]=max_nor_strength;
	(*shr_strength)[cmpnt]=max_shr_strength;
	(*ind_crack)[cmpnt] = 0;
      }
    }
    pos +=dx[0];
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Interfacer<_linear_coupled_cohesive>::createRightPropagatingCrackRoundAsp(Real initial_crack_size,
										      Real crit_nor_opening, 
										      Real crit_shr_opening,
										      Real max_nor_strength, 
										      Real max_shr_strength, 
										      Real radius, 
										      std::vector<Real> asp_ctr,
										      Real ratio) {

  std::vector<Real> pos = {0.5*dx[0], 0.5*dx[1]}; 

  crit_n_open.resize(total_n_ele, crit_nor_opening);
  crit_s_open.resize(total_n_ele, crit_shr_opening);

  for (UInt i = 0; i < n_ele[0]; ++i) {
    pos[1]=0.5*dx[1];
    for (UInt j = 0; j < n_ele[1]; ++j) {
     
      UInt cmpnt = i+n_ele[0]*j;

      if (pos[0] <= initial_crack_size){
	(*nor_strength)[cmpnt]=0;
	(*shr_strength)[cmpnt]=0;
	(*ind_crack)[cmpnt] = 2;
      }
           
      //round Asperity
      else if ((pos[0]-asp_ctr[0])*(pos[0]-asp_ctr[0])+(pos[1]-asp_ctr[1])*(pos[1]-asp_ctr[1])<radius*radius){

	(*nor_strength)[cmpnt]=ratio*max_nor_strength;
	(*shr_strength)[cmpnt]=ratio*max_shr_strength;
	
	(*ind_crack)[cmpnt]=5;
	}
      
      else{
	
	(*nor_strength)[cmpnt]=max_nor_strength;
	(*shr_strength)[cmpnt]=max_shr_strength;
	(*ind_crack)[cmpnt] = 0;
      }
   
      pos[1] += dx[1];
    }
    pos[0] += dx[0];
  }
}

/* -------------------------------------------------------------------------- */
template<>
void Interfacer<_linear_coupled_cohesive>::createIncohIntfc() {

  for (UInt i = 0; i < total_n_ele; ++i) {
    
      (*nor_strength)[i]=0;
      (*shr_strength)[i]=0;
      (*ind_crack)[i] = 2;
  }
}
