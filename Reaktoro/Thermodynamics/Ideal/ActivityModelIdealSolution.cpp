// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#include "ActivityModelIdealSolution.hpp"

namespace Reaktoro {

auto ActivityModelIdealSolution(StateOfMatter stateofmatter) -> ActivityModelGenerator
{
    ActivityModelGenerator model = [=](const SpeciesList& species)
    {
        ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args)
        {
            // Set the state of matter of the phase
            props.som = stateofmatter;

            props = 0.0;
            props.ln_a = args.x.log();
        };

        return fn;
    };

    return model;
}

} // namespace Reaktoro
