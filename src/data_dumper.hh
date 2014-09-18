/**
 * @file   data_dumper.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Sep  3 17:04:05 2014
 *
 * @brief  Class dealing with the outputs generation
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

/* -------------------------------------------------------------------------- */
#ifndef __DATA_DUMPER__
#define __DATA_DUMPER__
/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
#include "crack_profile.hh"
#include <iostream>
#include <map>
/* -------------------------------------------------------------------------- */
enum DumperType {
  _st_diagram,
  _snapshot,
  _points_history,
};

enum STDiagramType {
  _cracking_index,
  _shear_traction,
  _normal_traction
};

enum SnapshotType {
  _top_displacement,
  _bot_displacement,
  _tractions
};

/* -------------------------------------------------------------------------- */
class STDiagramBuilder {
public:
  STDiagramBuilder(const std::string & filename, int sze) {
    file.open(filename.c_str());
    size = sze;
  }

  virtual void dump(int stride, int start) = 0;

  std::ofstream file;

protected:
  int size;
};

/* -------------------------------------------------------------------------- */
template<class Bed>
class STDiagramDumper : public STDiagramBuilder {
public:
  STDiagramDumper(const std::string & filename, int nb_elem, const Bed & bed) :
    STDiagramBuilder(filename,nb_elem), data(bed) {}
  ~STDiagramDumper() {if (file.is_open()) file.close();}
  void dump(int stride, int start);

private:
  const Bed & data;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
class DataDumper {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  DataDumper(SpectralModel & mdl, const std::string sim_description, 
	     std::string release_info);

  virtual ~DataDumper();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // Initialize the energetic outputs
  void initEnergetics(const std::string filename);
  // Initialize the space time diagram output of a given type
  void initSpaceTimeDiagram(const std::string filename, STDiagramType type);
  // Initialize interface snap shot outputs
  void initSnapshot(const std::string filename, SnapshotType type);
  // Initialize point history outputs at given positons along the interface
  void initPointsHistory(const std::string filename, std::vector<double> positions);
  // Dump energetic data at time step it
  void printEnergetics(int it);
  // Dump at time step it every space-time diagrams defined
  void printSpaceTimeDiagram(int it);
  // Dump at time step it every snap shot plots defined
  void printSnapshot(int it);
  // Dump history of interface parameters at different positions of the interface
  void printPointsHistory(int it);

private:

  // Initialize displacements and velocities pointers
  void initDisplacementFields();
  // Initialize tractions and strengths pointers
  void initTractionFields();
  // Initialize the timer output related to the type of dump
  void initTimer(DumperType type);
  // Dump time step it
  void printTime(int it, DumperType type);
  // Construct summary files of simulation parameters
  void printSimulationSummary(const std::string sim_description);
  // Return the current time
  std::string getCurrentTime();

public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
 
  // Model generating the data dumped
  SpectralModel * model;
  // Info on the sources generated at compilation
  std::string cR_release_info;
  // Please refer to spectral_model.hh for info about these model members
  double beta;
  double dx;
  double X;
  int n_ele;
  int dim;
  const std::vector<Energetics> * E_nor;
  const std::vector<Energetics> * E_shr;
  const std::vector<Energetics> * E_fri;
  const std::vector<CrackProfile> * displacements;
  const std::vector<CrackProfile> * velocities;
  const CrackProfile * intfc_trac;
  const std::vector<double> * nor_strength;
  const std::vector<double> * shr_strength;
  const std::vector<unsigned int> * ind_crack;
  InterfaceFields * displ_jump;
  InterfaceFields * veloc_jump;
  // Array containing the point where history should be dumped
  std::vector<double> points;
  // Ofstream for summary file
  std::ofstream summary;
  // Ofstream for energetic outputs
  std::ofstream E_file;
  // Ofstream for point history outputs
  std::ofstream phistory_file;
  // Array containing the STDiagramDumper of the different ST diagrams
  std::map<STDiagramType,STDiagramBuilder*> st_diagram_dumper;
  // Array containing the STDiagramDumper of the different snap shot plots
  std::map<SnapshotType,STDiagramBuilder*> snapshot_dumper;
  // Map between a timer and its related dumper type
  std::map<DumperType,std::ofstream*> timer;
  // Map between ofstreams pointer and filename
  std::map<std::ofstream*, std::string> output_names ;
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "data_dumper_inline_impl.cc"

#endif /* __DATA_DUMPER__ */
