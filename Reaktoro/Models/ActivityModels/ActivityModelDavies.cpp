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

#include "ActivityModelDavies.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::log;
using std::pow;
using std::sqrt;

namespace detail {

/// Return the ActivityModel object based on the Davies model.
auto activityModelDavies(const SpeciesList& species, ActivityModelDaviesParams params) -> ActivityModel
{
    // Create the aqueous mixture
    AqueousMixture mixture(species);

    // The molar mass of water
    const auto Mw = mixture.water().molarMass();

    // The number of moles of water per kg
    const auto nwo = 1.0/Mw;

    // The number of charged and neutral species in the aqueous mixture
    const auto num_charged_species = mixture.charged().size();
    const auto num_neutral_species = mixture.neutral().size();

    // The indices of the charged and neutral species
    const auto icharged_species = mixture.indicesCharged();
    const auto ineutral_species = mixture.indicesNeutral();

    // The index of the water species
    const auto iwater = mixture.indexWater();

    // The electrical charges of the charged species only
    const ArrayXd charges = mixture.charges()(icharged_species);

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto stateptr = std::make_shared<AqueousMixtureState>();
    auto mixtureptr = std::make_shared<AqueousMixture>(mixture);

    // Define the activity model function of the aqueous mixture
    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the aqueous mixture
        auto const& state = *stateptr = mixture.state(T, P, x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;

        // Export the aqueous mixture and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = stateptr;
        props.extra["AqueousMixture"] = mixtureptr;

        // Auxiliary constant references
        const auto& m = state.m;             // the molalities of all species
        const auto& ms = state.ms;           // the stoichiometric molalities of the charged species
        const auto& I = state.Is;            // the stoichiometric ionic strength
        const auto& rho = state.rho/1000;    // the density of water (in g/cm3)
        const auto& epsilon = state.epsilon; // the dielectric constant of water

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Auxiliary variables
        const auto ln_m = m.log();
        const auto xw = x[iwater];
        const auto ln_xw = log(xw);
        const auto I2 = I*I;
        const auto sqrtI = sqrt(I);
        const auto sqrt_rho = sqrt(rho);
        const auto T_epsilon = T * epsilon;
        const auto sqrt_T_epsilon = sqrt(T_epsilon);
        const auto A = 1.824829238e+6 * sqrt_rho/(T_epsilon*sqrt_T_epsilon);
        const auto bions = params.bions;
        const auto bneutrals = params.bneutrals;
        const auto sigmac = -A*(sqrtI/(1 + sqrtI) - bions*I) * ln10;
        const auto sigman = bneutrals*I * ln10;

        // Calculate the contribution of ions to the ln activity of water
        ln_a[iwater] = ln10 * Mw * A * (2*(I + 2*sqrtI)/(1 + sqrtI) - 4*log(1 + sqrtI) - bions*I2) - (1 - xw)*xw;

        // Loop over all charged species in the aqueous mixture
        for(Index i = 0; i < num_charged_species; ++i)
        {
            // The index of the current charged species
            const auto ispecies = icharged_species[i];

            // The electrical charge of the charged species
            const auto zi = charges[i];

            // Calculate the ln activity coefficient of the current charged species
            ln_g[ispecies] = sigmac * zi*zi;

            // Calculate the ln activity of the current charged species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];
        }

        // Loop over all neutral species in the aqueous mixture
        for(Index i = 0; i < num_neutral_species; ++i)
        {
            // The index of the current neutral species
            const auto ispecies = ineutral_species[i];

            // The integral term Γᵢ for the neutral species
            const auto Gammai = m[ispecies] * sigman;

            // Calculate the ln activity coefficient of the current neutral species
            ln_g[ispecies] = sigman;

            // Calculate the ln activity coefficient of the current neutral species
            ln_a[ispecies] = ln_g[ispecies] + ln_m[ispecies];

            // Calculate the contribution of current neutral species to the ln activity of water
            ln_a[iwater] -= Mw * (m[ispecies]*ln_g[ispecies]);
        }

        // Set the activity coefficient of water (mole fraction scale)
        ln_g[iwater] = ln_a[iwater] - ln_xw;
    };

    return fn;
}

} // namespace detail

auto ActivityModelDavies() -> ActivityModelGenerator
{
    return ActivityModelDavies(ActivityModelDaviesParams{});
}

auto ActivityModelDavies(ActivityModelDaviesParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelDavies(species, params);
    };
}

} // namespace Reaktoro
