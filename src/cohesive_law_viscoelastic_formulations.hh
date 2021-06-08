/**
 * @file   cohesive_law_viscoelastic_formulations.hh
 * @author Thibault Roch <thibault.roch@epfl.ch>
 * @date   Wed Nov 18 18:56:06 2020
 *
 * @brief  Struct describing the formulations of the viscoelastic laws
 *
 * @section LICENSE
 *
 * cRacklet - A spectral boundary integral method for interface fracture simulation
 * Copyright (©) 2012 - 2013 Fabian Barras
 *               2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 *               Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 */
/* -------------------------------------------------------------------------- */
struct ViscoelasticFormulation {
  ViscoelasticFormulation(){};
  virtual ~ViscoelasticFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0*(1+std::max(0.,rate)/vel_lim);
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      slope = strength_v0/vel_lim;
    }
    return slope;
  }
};

/* -------------------------------------------------------------------------- */
struct ViscoelasticQuadraticFormulation : public ViscoelasticFormulation{
  ViscoelasticQuadraticFormulation(){};
  virtual ~ViscoelasticQuadraticFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    return strength_v0*(1+(std::max(0.,rate)/vel_lim)*(std::max(0.,rate)/vel_lim));
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      slope = 2*strength_v0*rate/(vel_lim*vel_lim);
    }
    return slope;
  }
};

/* -------------------------------------------------------------------------- */
struct ViscoelasticPowerLawFormulation : public ViscoelasticFormulation {
  ViscoelasticPowerLawFormulation(){};
  virtual ~ViscoelasticPowerLawFormulation(){};
  virtual inline Real getStrength(Real strength_v0, Real rate, Real vel_lim){
    Real exp = 1 / ( 1 - std::abs(rate)/vel_lim);
    return pow(strength_v0,exp);
  }
  virtual inline Real getTangent(Real strength_v0, Real rate, Real vel_lim){
    Real slope;
    if(rate < 0){
      slope = 0.;
    }else{
      Real exp = 1 / ( 1 - std::abs(rate)/vel_lim);
      slope = vel_lim * log(strength_v0) * pow(strength_v0,exp) / pow((vel_lim - std::abs(rate)),2);
    }
    return slope;
  }
};
