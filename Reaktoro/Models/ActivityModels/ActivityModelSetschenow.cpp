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

#include "ActivityModelSetschenow.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>

namespace Reaktoro {

using std::log;

auto ActivityModelSetschenow(String neutral, real b) -> ActivityModelGenerator
{
    ActivityModelGenerator model = [=](const SpeciesList& species)
    {
        // The index of the neutral aqueous species in the aqueous phase.
        const auto ineutral = species.indexWithFormula(neutral);

        ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args)
        {
            // Check AqueousMixtureState is available in props.extra
            auto stateit = props.extra.find("AqueousMixtureState");

            errorif(stateit == props.extra.end(),
                "ActivityModelSetschenow expects that another aqueous activity model has been chained first (e.g., Davies, Debye-Huckel, HKF, PitzerHMW, etc.) ");

            // The aqueous mixture state exported by a base aqueous activity model.
            const auto& state = *std::any_cast<SharedPtr<AqueousMixtureState> const&>(stateit->second);

            const auto& I = state.Is;
            props.ln_g[ineutral] = ln10 * b * I;
            props.ln_a[ineutral] = props.ln_g[ineutral] + log(state.m[ineutral]);
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
