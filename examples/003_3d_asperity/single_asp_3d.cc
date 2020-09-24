/**
 * @file   alu_homa3d.cc
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @date   Mon Nov 16 10:03:54 2015
 *
 * @brief  3d interface in the presence of one asperity (See http://lsms.epfl.ch/page-74218-en.html)
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation 
 * Copyright (©) 2012 - 2013 Fabian Barras 
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne) 
 *  
 * cRacklet is the result of a collaboration between the Computational Solid Mechanics  
 * Laboratory (LSMS) of Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland  
 * and the Department of Aerospace Engineering of the University of Illinois at  
 * Urbana-Champaign, United States of America. 
 *  
 * cRacklet is free software: you can redistribute it and/or modify it under the terms  
 * of the GNU General Public License as published by the Free Software Foundation,  
 * either version 3 of the License, or (at your option) any later version. 
 *  
 * cRacklet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;  
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   
 * See the GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License along with this program.   
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "simulation_driver.hh"
#include "interfacer.hh"
#include "cohesive_law.hh"
#include "coulomb_law.hh"
#include "regularized_coulomb_law.hh"
#include "data_dumper.hh"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <string>
#include <iomanip>
#include <sys/stat.h>
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]){

  // Note : Construct the pre-integrated material kernels before running this simulation
  // Use "invert_serial.f" to construct kernel files

  std::cout << "./single_asp_3d <output_folder_name> <asperity_toughness ratio> <nb_ele_x> <nb_time_steps> [ <asperity_x_center*G_length=2.25> <loading_angle=0.0> ]" << std::endl;

  std::time_t start = std::time(NULL);
  Real duration;
  Real pausing_duration = 72*60*60;
  bool first_half = true;
  
  std::string output_folder=argv[1];
   
  // Geometry description
  UInt nb_time_steps = std::atoi(argv[4]); 
  UInt nex = std::atoi(argv[3]); 
  UInt nez = nex/4;
  Real nu =  0.35;
  Real E = 5.3e9;
  Real cs = 1263;

  // Cut of the loaded material kernels
  UInt tcut = 100; 
  
  // Loading case
  Real load = 1e6;
  Real psi =0.0;
  if(argc > 6)
    psi = std::atof(argv[6]);
  Real phi = 0.0;

  // Cohesive parameters
  Real crit_n_open = 0.02e-3*sqrt(10.);
  Real crit_s_open = 0.02e-3*sqrt(10.);
  Real max_n_str = 5e6/sqrt(10.);
  Real max_s_str = 5e6/sqrt(10.);
  Real G_length = crit_n_open*max_n_str/(load*load*M_PI)*E/(1-nu*nu);

  Real dom_sizex = 10*G_length;
  Real dom_sizez = 2.5*G_length;
  Real dx = dom_sizex/(Real)(nex);
  Real dz = dom_sizez/(Real)(nez);

  Real crack_size = 2*dx;
   
  //ratio asperity strength/surrounding material strength
  Real ratio = 3.;
  ratio = std::atof(argv[2]);
  Real asper_radius=0.5*G_length;
  Real asperx = 2.25;
  Real asperz = 0.5*dom_sizez;
   
  if(argc > 5)
    asperx = std::atof(argv[5]);

  asperx*=G_length;

  // Friction parameters
  bool overlap = 0;
   
  std::string sim_name = "Mode-I fracture along a 3d interface with a single asperity";

  std::cout << "./single_asp_3d " 
	    << "output folder: " << output_folder << " " 
	    << "toughness ratio: " << ratio << " " 
	    << "nb_elements alog x: " << nex << " "
	    << "nb_time_steps: " << nb_time_steps << " "
	    << "asperity center position (x-coord/G_length): " << asperx/G_length << " "
	    << "griffith crack length: " << G_length << " "
	    << std::endl;

   
  /* -------------------------------------------------------------------------- */

  UInt t = 0;
  UInt x_tip=0;
  UInt x_lap = 0.05*nex;

  std::string outfolder;

  UInt nb_simulation_phases;
  if(ratio==1)
    nb_simulation_phases=1;
  else
    nb_simulation_phases=2;
  
  for (UInt phase = 0; phase < nb_simulation_phases; ++phase) {
    
    SpectralModel * model;
    
    if ((ratio==1.)||(phase==0)){
      outfolder = output_folder+"2d_outputs/";
      mkdir(outfolder.c_str(),0777);

#if defined (_OPENMP)
      max_num_threads = omp_get_max_threads();
      omp_set_num_threads(1);
#endif
      model = new SpectralModel(nex, 0, dom_sizex,
				nu, E, cs, tcut,
				sim_name, outfolder);      
    } else {
      outfolder = output_folder;
      mkdir(outfolder.c_str(),0777);

#if defined (_OPENMP)
      omp_set_num_threads(max_num_threads);
#endif
      model = new SpectralModel({nex,nez}, 0, {dom_sizex,dom_sizez},
				nu, E, cs, tcut,
				sim_name, outfolder);
    }
    
    SimulationDriver sim_driver(*model);
   
    Interfacer<_linear_coupled_cohesive> interfacer(*model);   

    if ((ratio==1.)||(phase==0)) {
      DataRegister::registerParameter("critical_normal_opening",crit_n_open);
      DataRegister::registerParameter("critical_shear_opening",crit_s_open);
      DataRegister::registerParameter("max_normal_strength",max_n_str);
      DataRegister::registerParameter("max_shear_strength",max_s_str);
      interfacer.createUniformInterface();
      interfacer.createThroughCrack(0.,crack_size);
    } else
      interfacer.createRightPropagatingCrackRoundAsp(crack_size, crit_n_open,crit_s_open,max_n_str,
						     max_s_str,asper_radius, {asperx,asperz},
						     sqrt(ratio),sqrt(ratio));

    interfacer.createThroughWall(0.95*dom_sizex,dom_sizex);
    CohesiveLaw& cohesive_law = dynamic_cast<CohesiveLaw&>((model->getInterfaceLaw()));
    cohesive_law.preventSurfaceOverlapping(NULL);
    
    sim_driver.initConstantLoading(load, psi, phi);
    
    /* -------------------------------------------------------------------------- */
    //Set-up simulation  outputs
     
    DataDumper dumper(*model);
     
    Real x_start = asperx-asper_radius; //0.4;
    Real x_end = x_start + 2.5*G_length; //1.0;
    UInt step = 5;
    std::vector<Real> z_coord;
     
    if ((ratio==1.)||(phase==0))
      z_coord = {0.};
    else
      z_coord = {1.25*G_length, 1.5*G_length, 1.75*G_length, 2.0*G_length}; //{0.25, 0.3, 0.35, 0.4};
    
    std::vector<UInt> obs_points;
    
    for (UInt z = 0; z < z_coord.size(); ++z) {
      UInt iz = (z_coord[z]/dz);
      UInt ix = (x_start/dx);
      while (dx*(Real)ix<x_end){
	obs_points.push_back(iz*nex+ix);
	ix += step;
      }
    }
    std::vector<DataFields> fields =
      {_interface_tractions,_normal_strength,_normal_displacement_jumps,_normal_velocity_jumps};
    
    std::cout << obs_points.size() << " observation points used to monitor fields evolution." << std::endl;   

    if(phase==1)
      dumper.initPointsDumper("mid-points_fields.cra", fields, obs_points, _binary);
    
    dumper.initDumper("ST_Diagram_nor_velo_jump.cra", _normal_velocity_jumps, 1.0, 1, 0, _binary);
    dumper.initDumper("ST_Diagram_id.cra", _id_crack, 1.0, 1, 0, _binary);
    //dumper.initDumper("ST_Diagram_nor_strength.cra", _normal_strength, 0.25, 4, 0, _binary);
    //dumper.initVectorDumper("ST_Diagram_nor_traction.cra", _interface_tractions, 1, 0.25, 4, 0, _binary);
    //Real ratio = 1/((UInt)(nez));
    // dumper.initDumper("ST_Diagram_center_id.cra", _id_crack, ratio, 1, nex*nez/2, _binary);

    /* -------------------------------------------------------------------------- */
    
    if(phase==1) {
      DataRegister::restart_dir = "2d_outputs/restart_files/";
      model->restartModel(true);
    }
        
    UInt chkptx = (asperx-asper_radius)/dx;
    UInt chkptz;
    
    if (phase==0)
      chkptz = 0;
    else
      chkptz = asperz/dz;
    
    const CrackProfile * interface_tractions = model->readData(_interface_tractions);
    
    if(phase==0) {
      sim_driver.launchCrack(0.,1.75*G_length,0.025);
    }
    
    while ((t < nb_time_steps)&&(x_tip<0.55*nex)) { //65

      if (phase==0) {
	model->pauseModel();
	std::cout << "End of pseudo 2d" << std::endl;
	break;
      }

      sim_driver.solveStep();
      x_tip = model->getCrackTipPosition(0.,nex);
      
      Real asperity_trac = (*interface_tractions)[(chkptx+nex*chkptz)*3+1];
      
      if (t%10==0){
	dumper.dumpAll();
      }
      else if (phase==1) {
	dumper.dump("mid-points_fields.cra");
	//dumper.dump("ST_Diagram_center_id.cra");
      }
      
      if ((x_tip>x_lap)||(t%(UInt)(0.05*nb_time_steps)==0)) {
	std::cout << "Process at " << (Real)t/(Real)nb_time_steps*100 << "% " << std::endl;
	std::cout << "Crack at " << 100*x_tip/(Real)(nex) << "% " << std::endl;
	std::cout << "Tractions at the tip of the asperity "
		  << asperity_trac << " [Pa]" << std::endl;
	std::cout << std::endl;
	
	if (x_tip>x_lap)
	  x_lap += 0.05*nex;
      }
      ++t;

      duration = std::time(NULL) - start;
      
      if((duration>0.5*pausing_duration)&&(first_half)) {
	first_half=false;
	std::cout << "The simulation just passes half of his planned duration !" <<std::endl;
      }
      
      /*      if ((phase==0)&&(asperity_trac>0.5*max_n_str)&&(ratio!=1.)) {
	      model->pauseModel();
	      std::cout << "End of pseudo 2d" << std::endl;
	      break;
	}
      */	
    }
    delete model;
  }
  return 0;
}

