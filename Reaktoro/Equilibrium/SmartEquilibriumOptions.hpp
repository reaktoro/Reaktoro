// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>

namespace Reaktoro {

/// The options for the smart equilibrium calculations.
/// @see SmartEquilibriumSolver
struct SmartEquilibriumOptions
{
    /// The options for the chemical equilibrium calculations during learning operations.
    EquilibriumOptions learning;

    /// The relative tolerance for negative species amounts when predicting with first-order Taylor approximation.
    /// When performing a first-order Taylor approximation for the calculation of a chemical
    /// equilibrium state, negative amounts of species may result. However, as long as these
    /// negative values are small enough, we should not rule out the approximate chemical
    /// equilibrium state. This relative tolerance value can be set to control the decision to
    /// accept or not the predicted chemical state with negative species amounts. If the accepted
    /// predicted equilibrium state has negative amounts of species, these are modified and set to
    /// their respective positive lower bounds.
    double reltol_negative_amounts = -1.0e-14;

    /// The relative tolerance used for the amount of species in the conservation of mass constraint.
    double reltol_component_amount_conservation = 1.0e-10;

    /// The relative tolerance used in the acceptance test for the predicted chemical equilibrium state.
    double reltol = 0.005;

    /// The absolute tolerance used in the acceptance test for the predicted chemical equilibrium state.
    double abstol = 0.01;

    /// The step length used to discretize temperature in the temperature-pressure space when storing learned calculations (in K).
    double temperature_step = 10.0;

    /// The step length used to discretize pressure in the temperature-pressure space when storing learned calculations (in Pa).
    double pressure_step = 25.0e+5;
};

} // namespace Reaktoro
