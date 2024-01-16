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

#include "ReactionStandardThermoModelPressureCorrection.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

auto ReactionStandardThermoModelPressureCorrection(real const& Pr) -> ReactionStandardThermoModel
{
    auto evalfn = [=](ReactionStandardThermoProps& props, ReactionStandardThermoModelArgs args)
    {
        // Unpack the arguments for the evaluation of this model
        const auto& [T, P, dV0] = args;

        const auto dE = (P - Pr) * dV0; // delta energy (in J/mol)
        props.dG0 += dE;
        props.dH0 += dE;
        props.dCp0 += 0.0; // TODO: Consider adding (P - Pr) * ΔVT0 in the future, taking into account temperature derivative of dV0.
    };

    Vec<real> params = { Pr };

    return ReactionStandardThermoModel(evalfn, params);
}

} // namespace Reaktoro
