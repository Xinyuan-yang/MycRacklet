/* -------------------------------------------------------------------------- */
template<class Bed>
void STDiagramDumper<Bed>::dump(int stride, int start) {

  for (int i = 0; i < size; ++i) {

    file << std::scientific << data[i*stride+start] << " "; 
 }
  file << std::endl;
}

/* -------------------------------------------------------------------------- */
inline DataDumper::DataDumper(SpectralModel & mdl, const std::string sim_description,
			      std::string release_info) {

  model = & mdl;
  
  std::vector<double> sim_parameters = model->getSimulationsParameters();
  X = sim_parameters[0];
  dx = X/sim_parameters[1];
  beta = sim_parameters[2]; 
  n_ele = (int)sim_parameters[1];
  dim = (int)sim_parameters[3];

  cR_release_info = release_info;

  E_nor = NULL;
  E_shr = NULL;
  E_fri = NULL;
  displacements = NULL;
  velocities = NULL;
  intfc_trac = NULL;
  nor_strength = NULL;
  shr_strength = NULL;
  ind_crack = NULL; 

  displ_jump = NULL;
  veloc_jump = NULL;

  printSimulationSummary(sim_description);

}
/* -------------------------------------------------------------------------- */
inline DataDumper::~DataDumper(){

  if(phistory_file.is_open()) phistory_file.close();

  if(displ_jump) delete displ_jump;
  if(veloc_jump) delete veloc_jump;

  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " GENERATED OUTPUT FILES " << std::endl;

  std::string filename;

  std::map<DumperType,std::ofstream*>::iterator it = timer.find(_st_diagram);

  if(E_file.is_open()) {
   
    summary << "* Energetics data: " << std::endl 
	    << output_names[&E_file] << std::endl;
    E_file.close();
  }
  
  if (it != timer.end()) {
    summary << std::endl;
    summary << "* Space Time Diagram:  (Timed by " << output_names[timer[_st_diagram]]  
	    << " )" << std::endl;
    
    (it->second)->close();
    delete it->second;
    timer.erase(it);
  }
  
  std::map<STDiagramType,STDiagramBuilder*>::iterator itt;
  for (itt = st_diagram_dumper.begin(); itt != st_diagram_dumper.end(); ++itt){  

    filename = output_names[&((itt->second)->file)];
    summary << filename;
    summary << std::endl;
    delete itt->second;

  }

  it = timer.find(_snapshot);
  if (it != timer.end()) {
    summary << std::endl;
    summary << "* Interface Snapshot:  (Timed by " << output_names[timer[_snapshot]]  
	    << " )" << std::endl;
    
    (it->second)->close();
    delete it->second;
    timer.erase(it);
  }

  std::map<SnapshotType,STDiagramBuilder*>::iterator ittt;
  for (ittt = snapshot_dumper.begin(); ittt != snapshot_dumper.end(); ++ittt){  

    filename = output_names[&((ittt->second)->file)];
    summary << filename;
    summary << std::endl;
    delete ittt->second;
  }

  it = timer.find(_points_history);
  if (it != timer.end()) {
    summary << std::endl;
    summary << "* Failure history at given positions:  (Timed by " << output_names[timer[_points_history]] 
	     << " )" << std::endl;
    (it->second)->close();
    delete it->second;
    timer.erase(it);
  }

  filename = output_names[&phistory_file];
  summary << filename;
  summary << std::endl;


  for (it = timer.begin(); it != timer.end(); ++it){  
    
    (it->second)->close();
    delete it->second;

  }

  summary << std::endl;
  summary << "/* -------------------------------------------------------------------------- */ "; 
  summary << std::endl;
  summary << " Thank you for using cRacklet, hope to see you soon.";
 
  summary.close();

}
