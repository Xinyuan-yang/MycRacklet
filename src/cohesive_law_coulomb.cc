/**
 * @file   cohesive_law_coulomb.cc
 * @author Roxane Ferry <roxane.ferry@epfl.ch>
 * @date   Tue Oct 3 11:16:47 2023
 *
 * @brief  Implementation of the CohesiveLawCoulomb class
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
#include "cohesive_law_coulomb.hh"
#include <fstream>
#include <algorithm>
#include <random>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::initStandardFormulation()
{
  formulation = std::make_shared<CohesiveCoulombFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::initDualFormulation()
{
  formulation = std::make_shared<DualCohesiveCoulombFormulation>();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::initInterfaceConditions()
{
  // memcpy(&this->cf[0],&this->cf_s[0],this->cf_s.size()*sizeof(Real));
  /*
  UInt n_ele = cf_s.size();
  for (UInt i = 0; i < n_ele; ++i) {
    fric_strength[i] = sigma_0*cf_s[i];
  }
  */
  this->computeInitialVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::updateInterfaceConditions()
{
  this->updateCohesiveLaw();
  this->computeVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::updateCohesiveLaw()
{
  CrackProfile *shr_opening = datas[_shear_displacement_jumps];
  CrackProfile *nor_opening = datas[_normal_displacement_jumps];
  CrackProfile *loads = datas[_top_loading];

  UInt n_ele = cf_s.size();
  Real aux;

  for (UInt i = 0; i < n_ele; ++i)
  {
    aux = sqrt(((*nor_opening)[i] / crit_nor_opening[i]) * ((*nor_opening)[i] / crit_nor_opening[i]) +
               ((*shr_opening)[i] / crit_shr_opening[i]) * ((*shr_opening)[i] / crit_shr_opening[i]));
    if (((*nor_opening)[i] == 0) && ((*shr_opening)[i] == 0) && (ind_crack[i] != 2))
    {
      cf[i] = cf[i]; // cf_s[i];
    }

    else if ((aux >= 1) || (fric_strength[i] == -(*loads)[i * 3 + 1] * cf_d[i])|| ind_crack[i] == 2 )
    {
      ind_crack[i] = 2;
      cf[i] = cf_d[i];
    }

    else
    {
      // ind_crack[i] = 1;
      cf[i] = (*formulation)(cf_s[i], cf_d[i], aux, (*shr_opening)[i], i, ind_crack[i]);
    }
    fric_strength[i] = std::max(-(*loads)[i * 3 + 1] * cf[i],0.0);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::restart(bool pausing, UInt nele_2d)
{
  std::vector<Real> *fric_strength = datas[_frictional_strength];
  DataRegister::restartData(*fric_strength, "restart_frictional_strength.cra", pausing, nele_2d);
  std::vector<Real> *cf_s = datas[_static_friction_coefficient];
  DataRegister::restartData(*cf_s, "restart_static_friction_coefficient.cra", pausing, nele_2d);
  std::vector<Real> *cf_d = datas[_dynamic_friction_coefficient];
  DataRegister::restartData(*cf_d, "restart_dynamic_friction_coefficient.cra", pausing, nele_2d);
  // For DualCohesiveCoulombFormulation
  if (DataRegister::hasParameter("int_friction_coefficent"))
  {
    std::vector<Real> *cf_int = datas[_int_friction_coefficient];
    DataRegister::restartData(*cf_int, "restart_int_friction_coefficient.cra", pausing, nele_2d);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::computeInitialVelocities()
{
  // Get the stresses
  std::vector<CrackProfile *> loads = {datas[_top_loading], datas[_bottom_loading]};

  for (UInt side = 0; side < 2; ++side)
  {
    for (UInt i = 0; i < n_ele[0]; ++i)
    {
      for (UInt j = 0; j < n_ele[1]; ++j)
      {

        (*stresses[side])[(i * dim + 1) + j * n_ele[0] * dim] = (*loads[side])[(i * dim + 1) + j * n_ele[0] * dim];
      }
    }
  }

  this->computeVelocities();
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::computeVelocities()
{
  CrackProfile deltaStresses(n_ele, dim);
  std::vector<Real> temp_veloc(dim);
  Real trac;

  deltaStresses = (*stresses[0]) - (*stresses[1]);

  Real cste = 1 / (mu[0] / cs[0] + mu[1] / cs[1]);

  (*velocities[0]) = deltaStresses * cste;

  for (UInt i = 0; i < n_ele[0]; ++i)
  {
    for (UInt j = 0; j < n_ele[1]; ++j)
    {
      (*velocities[0])[(i * dim + 1) + j * n_ele[0] * dim] *= 1 / (mu[0] * eta[0] / cs[0] + mu[1] * eta[1] / cs[1]) * (1 / cste);
    }
  }
  (*velocities[1]) = (*velocities[0]);

  for (UInt i = 0; i < n_ele[0]; ++i)
  {
    for (UInt j = 0; j < n_ele[1]; ++j)
    {

      Real trac = (*stresses[0])[(i * dim + 1) + j * n_ele[0] * dim] - mu[0] * eta[0] * (*velocities[0])[(i * dim + 1) + j * n_ele[0] * dim] / cs[0];
      (*intfc_trac)[(i * dim + 1) + j * n_ele[0] * dim] = trac;
      computeShearVelocities(fric_strength[i + n_ele[0] * j], i + j * n_ele[0]);
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::computeShearVelocities(Real strength, UInt i)
{
  std::vector<Real> trac(2);
  Real shr_trac;

  for (UInt j = 0; j < 2; ++j)
  {

    trac[j] = (*stresses[0])[i * dim + 2 * j] - mu[0] * (*velocities[0])[i * dim + 2 * j] / cs[0];
  }

  shr_trac = sqrt((trac[0] * trac[0]) + (trac[1] * trac[1]));

  if ((strength < shr_trac) || (strength == 0))
    computeIndepShearVelocities(strength, i);
  else
  {

    for (UInt j = 0; j < (dim - 1); ++j)
    {

      (*intfc_trac)[i * dim + 2 * j] = trac[j];
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveLawCoulomb::computeIndepShearVelocities(Real strength, UInt i)
{
  std::vector<Real> cmpted_stress(2);
  Real dyn_stress;
  Real shr_veloc;

  for (UInt side = 0; side < 2; ++side)
  {

    for (UInt j = 0; j < 2; ++j)
    {

      cmpted_stress[j] = (*stresses[side])[i * dim + 2 * j];
    }

    dyn_stress = sqrt((cmpted_stress[0] * cmpted_stress[0]) + (cmpted_stress[1] * cmpted_stress[1]));

    if (side == 0)
      shr_veloc = cs[0] / mu[0] * (dyn_stress - strength);

    else
      shr_veloc = cs[1] / mu[1] * (strength - dyn_stress);

    for (UInt j = 0; j < 2; ++j)
    {

      if (dyn_stress == 0)
      {
        (*velocities[side])[i * dim + 2 * j] = 0;
      }
      else
      {
        (*velocities[side])[i * dim + 2 * j] = shr_veloc * cmpted_stress[j] / dyn_stress;
      }
      if (side == 0)
      {
        if (dyn_stress == 0)
        {
          (*intfc_trac)[i * dim + 2 * j] = 0;
        }
        else
        {
          (*intfc_trac)[i * dim + 2 * j] = strength * cmpted_stress[j] / dyn_stress;
        }
      }
    }
  }
}