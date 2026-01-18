/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "config/config.hpp"

#include "Particle.hpp"
#include "random.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

/** Langevin thermostat for particle retarded force.
 *  @param[in]     jeffreys_langevin       Parameters
 *  @param[in]     p              Particle
 *  @param[in]     time_step      Time step
 *  @param[in]     kT             Thermal energy
 */
inline Utils::Vector3d retarded_friction_thermo_langevin(
    JeffreysLangevinThermostat const &jeffryes_langevin, Particle const &p,
    double time_step, double kT) {
  using namespace Thermostat;

  auto const pref_friction = jeffryes_langevin.pref_retarded_friction;
  auto const pref_noise = jeffryes_langevin.pref_noise_retarded;

  auto const friction_op = handle_particle_anisotropy(p, pref_friction);
  auto const noise_op = handle_particle_anisotropy(p, pref_noise);
  return friction_op * p.v() +
         noise_op * Random::noise_uniform<RNGSalt::JEFFREYS_LANGEVIN>(
                        jeffryes_langevin.rng_counter(),
                        jeffryes_langevin.rng_seed(), p.id());
}

inline void retarded_force_propogator_half_step(
    JeffreysLangevinThermostat const &jeffreys_langevin, Particle &p,
    double time_step, double kT) {
  auto const retarded_friction_thermo =
      retarded_friction_thermo_langevin(jeffreys_langevin, p, time_step, kT);
  auto const relax_time = jeffreys_langevin.relax_time;

  p.retarded_f() += 0.5 * time_step *
                    -(p.retarded_f() - retarded_friction_thermo) / relax_time;
}
