/**
 * @file   data_dumper.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Wed Sep  3 17:04:05 2014
 *
 * @brief  Classes dealing with outputs production from SpectralModel fields
 * @brief  Dumper are object used to produce the desired output while DataDumper manage them     
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
#include "cRacklet_common.hh"
#include "data_register.hh"
#include "spectral_model.hh"
#include "crack_profile.hh"
#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <vector>
#include <map>

/* -------------------------------------------------------------------------- */
// Format of generated output files: text or binary files
enum OutputFormat {
  _text,
  _binary
};

// A list of fields preset by default when using PointDumper
static std::vector<DataFields> standard_fields_history =
  {_normal_strength,_shear_strength,_interface_tractions,_normal_displacement_jumps,
   _shear_displacement_jumps, _normal_velocity_jumps, _shear_velocity_jumps,
   _top_velocities};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
// Dumper class managing simulation data output to file
class Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  Dumper(){};
  // Dumper is created providing field to dump, which amout (ratio),
  // the stride to appy while dumping and the first position to dump. 
  Dumper(DataFields field, Real data_ratio,
	 UInt data_stride, UInt data_start) {
    field_name=field;
    ratio = data_ratio;
    stride=data_stride;
    start=data_start;
  }
  ~Dumper() {if (file.is_open()) file.close();}
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */  
  // Initialize Dumper with the name of output file and its format
  void init(const std::string & filename, OutputFormat out_for=_text) {
    output_format = out_for;
    switch (out_for) {
    case _text:
      file.open(filename.c_str());
      break;
    case _binary:
      file.open(filename.c_str(),std::ios::out|std::ios::binary);
      break;
    default:
      cRacklet::error("DataDumper don't know how to open file using this OutputFormat");
      break;
    }
  }
  // Print a backspace in the output file
  void endl(){if(output_format==_text){file << std::endl;}}
  // Accessor to Dumper's ofstream
  std::ofstream & getFile(){return file;}
  // Output current state of registered data
  virtual void dump();

protected:
  // Download the data of type Bed from DataRegister before dumping them
  template<class Bed>
  inline const Bed & getData();
  // Output current state of registered data of type Bed
  template<class Bed>
  inline void dump();

private:
  // Dumping methods for OutputFormat _text
  template<class Bed>
  inline void dump_text(Bed & data, UInt size);
  // Dumping methods for OutputFormat _binary
  template<class Bed>
  inline void dump_binary(Bed & data, UInt size);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  //Format use to output data in a file
  OutputFormat output_format;
  //Name of the output file
  std::ofstream file;
  //Name of the field dumped by the Dumper (use to access data through DataRegister) 
  DataFields field_name;
  //Ratio of the total number of data that will be dumped
  Real ratio;
  //First position to be dumped to the output file
  UInt start;
  //Stride to move within data while dumping
  UInt stride;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
// Daughter Dumper use to dump several fields at some specific interface positions
class PointsDumper: public Dumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  //Construct an object by specifying the fields and the position (btw 0 and total_n_ele)
  //where to output fields value 
  PointsDumper(std::vector<DataFields> & data_fields,
	       std::vector<UInt> & points_to_dump, UInt total_nb_elem) : 
    fields(data_fields), points(points_to_dump) {
    
    this->ratio = 1.0/(Real)(total_nb_elem);
    this->stride=1;
  }
  ~PointsDumper(){};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  // Output current state of registered data
  void dump();
  // Only text output is currently available for PointsDumper
  template<class Bed>
  inline void dump_text(Bed & data, UInt size);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  //Position of interest where to dump fields
  std::vector<UInt> & points;
  //Fields of interest to dump
  std::vector<DataFields> & fields;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
// Dumpers manager interfacing with users
class DataDumper {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  DataDumper(SpectralModel & mdl);
  virtual ~DataDumper();
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // Initialize the energetic outputs
  void initEnergetics(const std::string filename);
  // Create a Dumper to output a field of given type into filename
  // Eventually specify which amout (ratio) of the full data,
  // the stride to appy while dumping and the first position to dump.
  void initDumper(const std::string filename, DataFields type, Real ratio_of_nele=1.0,
		  UInt stride=1, UInt start=0, OutputFormat format=_text);
  // Create a Dumper to output vectorial field (of size n_ele*dim).
  // Options are the same than initDumper + one entry to specify which direction to dump  
  void initVectorDumper(const std::string filename, DataFields type, UInt dimension_to_dump,
			Real ratio_of_nele=1.0, UInt stride=1, UInt start=0, OutputFormat format=_text);
  //Create a PointsDumper to output a list of fields at precise interface position specified
  //in points_to_dump
  void initPointsDumper(const std::string filename, std::vector<DataFields> & fields,
			std::vector<UInt> & points_to_dump);
  //Create a PointDumper as above but using the standard output fields listed on top of this file
  void initPointsDumper(const std::string filename, std::vector<UInt> & points_to_dump);
  //Create a PointDumper of standard fields along nb_obs_points between elements start and end  
  void initPointsDumper(const std::string filename, UInt start, UInt end, UInt nb_obs_points);
  //Generate an output of current model state to given file 
  void dump(std::string filename);
  //Generate output from every created dumper at once
  void dumpAll();
  // Dump energetic data
  void printEnergetics();

private:

  //Initialize output timer related to a precise Dumper identified by its related filename
  void initTimer(std::string filenname);
  //Print current simulation time 
  void printTime(std::string filename);

public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
 
  //Model generating the data dumped
  SpectralModel * model;
  //Map containing all the created Dumpers identified by the filename they are writing
  std::map<std::string,Dumper*> dumper;
  // Map between the timer of a given output file and the ofstream use to dump time information
  std::map<std::string,std::ofstream*> timer;
  // Pointer to energetic quantity
  const std::vector<Energetics> * E_nor;
  const std::vector<Energetics> * E_shr;
  const std::vector<Energetics> * E_fri;
  const Energetics * Eq;
  // Ofstream for energetic outputs
  std::ofstream E_file;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "data_dumper_inline_impl.cc"

#endif /* __DATA_DUMPER__ */
