/* -------------------------------------------------------------------------- */
#include "data_dumper.hh"
#include <time.h>
/* -------------------------------------------------------------------------- */
void DataDumper::initEnergetics(const std::string filename) {

  E_nor = & model->getNormalDissipatedEnergy();
  E_shr = & model->getShearDissipatedEnergy();
  E_fri = & model->getFrictionalEnergy();

  E_file.open(filename.c_str());
  output_names[&E_file] = filename;
  std::cout << "Energetic values are stored in file " << filename << std::endl;  
}

/* -------------------------------------------------------------------------- */
void DataDumper::initSpaceTimeDiagram(const std::string filename, STDiagramType type) {

  switch (type) {

  case _cracking_index:
    if (!ind_crack) ind_crack = & model->getCrackingIndex();
    st_diagram_dumper[type] = new STDiagramDumper<std::vector<unsigned int> >(filename, n_ele, *ind_crack);
    output_names[&((st_diagram_dumper[type])->file)] = filename;
    std::cout << "Space-Time diagram of the cracking index is built in file " << filename << std::endl;  
    break;

  case _normal_traction:
  case _shear_traction:
    if (!intfc_trac) initTractionFields();
    st_diagram_dumper[type] = new STDiagramDumper<CrackProfile>(filename, n_ele, *intfc_trac);
    output_names[&((st_diagram_dumper[type])->file)] = filename;
    std::cout << "Space-Time diagram of the traction(" << type-1 << ") is built in file " 
	      << filename << std::endl;   
    break;
    
  default:
    std::cout << "*** The STDiagramType " << type 
	      << " is not defined. See data_dumper.hh" << std::endl;
    break;
  }

  std::map<DumperType,std::ofstream*>::iterator it = timer.find(_st_diagram);

  if (it == timer.end()) initTimer(_st_diagram);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initSnapshot(const std::string filename, SnapshotType type) {


  switch (type) {

  case _top_displacement:
    if(!displacements) initDisplacementFields();
    snapshot_dumper[type] = new STDiagramDumper<CrackProfile>(filename, dim*n_ele, (*displacements)[0]);
    output_names[&((snapshot_dumper[type])->file)] = filename;
    std::cout << "Snapshots of the top displacement fields is built in file "<< filename << std::endl;  
    break;

  case _bot_displacement:
    if(!displacements) initDisplacementFields();
    snapshot_dumper[type] = new STDiagramDumper<CrackProfile>(filename, dim*n_ele, (*displacements)[1]);
    output_names[&((snapshot_dumper[type])->file)] = filename;
    std::cout << "Snapshots of the bottom displacement fields is built in file "<< filename << std::endl;
    break;

  case _tractions:
    if(!intfc_trac) initTractionFields();
    snapshot_dumper[type] = new STDiagramDumper<CrackProfile>(filename, dim*n_ele, *intfc_trac);
    output_names[&((snapshot_dumper[type])->file)] = filename;
    std::cout << "Snapshots of the interface tractions is built in file "<< filename << std::endl;
    break;

  default:
    std::cout << "*** The SnapShotType " << type 
	      << " is not defined. See data_dumper.hh" << std::endl;
    break;
  }

  std::map<DumperType,std::ofstream*>::iterator it = timer.find(_snapshot);
  
  if (it == timer.end()) initTimer(_snapshot);
}

/* -------------------------------------------------------------------------- */
void DataDumper::initPointsHistory(const std::string filename, std::vector<double> positions) {
  
  if(!displacements) initDisplacementFields();
  if(!intfc_trac) initTractionFields();

  displ_jump = new InterfaceFields(displacements, dim);
  veloc_jump = new InterfaceFields(velocities, dim);

  points = positions;

  phistory_file.open(filename.c_str());

  std::map<DumperType,std::ofstream*>::iterator it = timer.find(_points_history);
  if (it == timer.end()) initTimer(_points_history);

  output_names[&phistory_file] = filename;

  std::cout << "Fracture history at different points along the interface is stored in "
	    << filename << std::endl;

}

/* -------------------------------------------------------------------------- */
void DataDumper::initTimer(DumperType type) {

  timer[type] = new std::ofstream;

  switch (type) {

  case _st_diagram:
    timer[type]->open("Time_stdiag.dat");
    output_names[timer[type]] = "Time_stdiag.dat" ;
    break;

  case _snapshot:
    timer[type]->open("Time_snapshot.dat");
     output_names[timer[type]] = "Time_snapshot.dat" ;
    break;

  case _points_history:
    timer[type]->open("Time_history.dat");
     output_names[timer[type]] = "Time_history.dat" ;
    break;
  }
}

/* -------------------------------------------------------------------------- */
void DataDumper::initDisplacementFields() {

  displacements = & model->getDisplacements();
  velocities = & model->getVelocities();

}

/* -------------------------------------------------------------------------- */
void DataDumper::initTractionFields(){

  intfc_trac = & model->getInterfaceTractions();
  nor_strength = & model->getNormalStrength();
  shr_strength = & model->getShearStrength();

}

/* -------------------------------------------------------------------------- */
void DataDumper::printEnergetics(int it) {

  E_file << std::scientific 
	 << it*beta*dx/X <<" "<< (*E_shr)[0].E+(*E_nor)[0].E <<" " 
	 <<(*E_shr)[0].E_dot_old+(*E_nor)[0].E_dot_old <<" "<< (*E_fri)[0].E <<" "
	 << (*E_fri)[0].E_dot_old <<" "<< (*E_shr)[1].E+(*E_nor)[1].E <<" "
	 << (*E_shr)[1].E_dot_old+(*E_nor)[1].E_dot_old <<" "<< (*E_fri)[1].E <<" "
	 << (*E_fri)[1].E_dot_old <<" "<< (*E_nor)[0].E <<" "<< (*E_nor)[0].E_dot_old <<" "
	 << (*E_shr)[0].E <<" "<< (*E_shr)[0].E_dot_old <<" "<<	(*E_nor)[1].E <<" "
	 << (*E_nor)[1].E_dot_old <<" "<< (*E_shr)[1].E <<" "<< (*E_shr)[1].E_dot_old
	 << std::endl;

}

/* -------------------------------------------------------------------------- */
void DataDumper::printSpaceTimeDiagram(int t) {

  std::map<STDiagramType,STDiagramBuilder*>::iterator it;
  for (it = st_diagram_dumper.begin(); it != st_diagram_dumper.end(); ++it){  
    
    switch (it->first) {

    case _cracking_index:
      (it->second)->dump(1,0);

      break;
      
    case _shear_traction:
      (it->second)->dump(dim,0);
      break;
      
    case _normal_traction:
      (it->second)->dump(dim,1);
      break;

    default:
      std::cout << "The STDiagramType " << it->first 
		<< " is not defined. See data_dumper.hh" << std::endl;
      break;
    }
  }

  printTime(t, _st_diagram);
}

/* -------------------------------------------------------------------------- */
void DataDumper::printSnapshot(int t) {

  std::map<SnapshotType,STDiagramBuilder*>::iterator it;
  for (it = snapshot_dumper.begin(); it != snapshot_dumper.end(); ++it){

      (it->second)->dump(1,0);
  }

  printTime(t, _snapshot);
}

/* -------------------------------------------------------------------------- */
void DataDumper::printPointsHistory(int it) {

  int nb_points = points.size();

  displ_jump->computeJumpFields();
  veloc_jump->computeJumpFields();

  for (int i = 0; i < nb_points; ++i) {

    int select = points[i]*n_ele;
    
    phistory_file << std::scientific 
		  << ((double)select+0.5)*dx << " " << (*nor_strength)[select] << " " 
		  << (*shr_strength)[select] << " " << (*intfc_trac)[select*dim] << " " 
		  << (*intfc_trac)[select*dim+1] << " " << (*intfc_trac)[select*dim+2] 
		  << " " << (*displ_jump).delta_fields[1][select] 
		  << " " << (*displ_jump).delta_fields[0][select] 
		  << " " << (*veloc_jump).delta_fields[1][select] 
		  << " " << (*veloc_jump).delta_fields[0][select] << " " ; 
  }
  
  phistory_file << std::endl;
  printTime(it, _points_history);
} 
/* -------------------------------------------------------------------------- */
void DataDumper::printTime(int it, DumperType type) {
  
  *(timer[type]) << std::scientific << it*beta*dx/X << std::endl;

}
/* -------------------------------------------------------------------------- */
void DataDumper::printSimulationSummary(const std::string sim_description) {

  std::ofstream parameters;

  parameters.open("Parameters.dat");

  summary.open("Simulation_Summary.ibs");

  summary << " ** This file is automatically generated by cRacklet DataDumper ** " << std::endl;
  summary << std::endl;
  summary << " Thank you for using cRacklet !" << std::endl;
  summary << " Hereunder is a summary of your last simulation " << std::endl;
  summary << std::endl;
  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " * Simulation Description: " << sim_description << std::endl;
  summary << std::endl;
  summary << " * Date and Time:  " << getCurrentTime() << std::endl;
  summary << " * Sources information: " << std::endl;
  summary << cR_release_info << std::endl;
 
  model->printSelf(parameters, summary);

  parameters.close();
}

/* -------------------------------------------------------------------------- */
std::string DataDumper::getCurrentTime() {
  
  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);

}


