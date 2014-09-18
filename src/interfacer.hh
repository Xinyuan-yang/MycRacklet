/**
 * @file   interfacer.hh
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @date   Thu Sep 11 13:17:04 2014
 *
 * @brief This class is used to set interface properties  
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
#ifndef __INTERFACER__
#define __INTERFACER____
/* -------------------------------------------------------------------------- */
#include "spectral_model.hh"
/* -------------------------------------------------------------------------- */
enum InterfaceType {
  _centered_crack,
  _left_sided_crack,
  _incoherent
};

class Interfacer {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  
  Interfacer(SpectralModel & model):
  shr_strength(model.getShearStrength()), 
  nor_strength(model.getNormalStrength()), 
  ind_crack(model.getCrackingIndex()) {
    n_ele = shr_strength.size();
    const std::vector<double> param = model.getSimulationsParameters();
    crack_size = param[4];
    dx = param[0]/param[1];
  }

  virtual ~Interfacer(){};
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  // create material for a centered crack
  void createCenteredCrack(double max_nor_strength, double max_shr_strength);
  // create material for a left-sided crack
  void createLeftSidedCrack(double max_nor_strength, double max_shr_strength);
  // create conditions of dynamic sliding after right impacting
  void createIncohIntfc();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  // Shear and normal strength arrays of the interface
  std::vector<double> & shr_strength;
  std::vector<double> & nor_strength;
  // Cracking index
  std::vector<unsigned int> & ind_crack;
  // Number of elements at the interface
  int n_ele;
  // Space step
  double dx;
  // Initial crack size over the domain size
  double crack_size;
  
};


#endif /* __INTERFACER______ */
