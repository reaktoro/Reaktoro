// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

#include "ActivityModelSurfaceComplexation.hpp"

// Reaktoro includes
#include <Reaktoro/Singletons/Elements.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousMixture.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Common/Constants.hpp>

namespace Reaktoro {

using std::sqrt;
using std::log;

// Auxiliary constants
const auto F = faradayConstant;
const auto R = universalGasConstant;

namespace detail {

/// Return the SurfaceComplexationActivityModel object assuming no electrostatic effects.
auto activityModelSurfaceComplexationNoDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The number of surface complexation species in the current phase
    const auto num_species = surface.species().size();

    // The charges of the surface complexation species
    ArrayXd z = surface.charges();

    // The state of the complexation surface
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        surface_state = surface.state(T, P, x);

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        ln_a = x.log();

        // Initialized the ln of activity coefficients of the surface complexation species
        ln_g = ArrayXr::Zero(num_species);

        // Calculate Davies only if the AqueousPhase has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Auxiliary constant references properties and variables
            const auto I = aqstate.Is;          // the stoichiometric ionic strength
            const auto sqrtI = sqrt(I);
            const auto ln10 = log(10);
            const auto Agamma = 0.5095;         // the Debye-Huckel parameter

            // Calculate the ln activity coefficient of the surface complexation species using the Davies activity model
            ln_g = ln10*(-Agamma*z*z*sqrtI/(1 + sqrtI) - 0.3*I);
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

/// Return the SurfaceComplexationActivityModel object assuming the presence the Diffuse Double Layer (DDL) model.
auto activityModelSurfaceComplexationWithDDL(const SpeciesList& species, ActivityModelSurfaceComplexationParams params) -> ActivityModel
{
    // Create the complexation surface
    ComplexationSurface surface = params.surface;

    // The number of surface complexation species in the current phase
    const auto num_species = surface.species().size();

    // The charges of the surface complexation species
    ArrayXd z = surface.charges();

    // The state of the complexation surface
    ComplexationSurfaceState surface_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the state of the surface complexation
        surface_state = surface.state(T, P, x);

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        // Calculate ln of activities of surfaces species as the ln of molar fractions
        ln_a = x.log();

        // Initialized the ln of activity coefficients of the surface complexation species
        ln_g = ArrayXr::Zero(num_species);

        // Auxiliary constant references properties
        real I;

        // Otherwise, calculate the stoichiometric ionic strength if the Aqueous State has been already evaluated
        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export aqueous mixture state via `extra` data member
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Fetch the stoichiometric ionic strength
            I = aqstate.Is;
        }

        // Auxiliary variables
        const auto sigma = surface_state.sigma;
        surface_state.updatePotential(I);
        const auto psi = surface_state.psi;

        // Export the surface complexation and its state via the `extra` data member
        props.extra["ComplexationSurfaceState"] = surface_state;
        props.extra["ComplexationSurface"] = surface;

        // Calculate ln of gamma according to the coulombic correction, Appelo etal (2005), (7.44), p. 334
        ln_g = z*F*psi/(R*T);

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };

    return fn;
}

/// Return the Diffuse Double Layer (DDL) activity model based on the Dzombak and Morel (1990) model
auto activityModelDDL(const SpeciesList& species, ActivityModelDDLParams params) -> ActivityModel
{
    // Create the aqueous ddl_mixture
    AqueousMixture ddl_mixture(species);

    // The number of all species and surface complexation species in the current exchange phase only
    const auto num_species = species.size();

    // The ddl_state of the aqueous ddl_mixture
    AqueousMixtureState ddl_state;

    // Define the activity model function of the surface complexation phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Evaluate the ddl_state of the aqueous ddl_mixture
        ddl_state = ddl_mixture.state(T, P, x);

        // Export the surface complexation and its ddl_state via the `extra` data member
        props.extra["DiffusiveLayerState"] = ddl_state;

        // Auxiliary references
        auto& ln_g = props.ln_g;
        auto& ln_a = props.ln_a;

        real I = 0;
        ArrayXr z_aq = ArrayXr::Zero(num_species);
        real psi;

        if (props.extra["AqueousMixtureState"].has_value())
        {
            // Export surface complexation ddl_state via `extra` data member
            const auto& aq_state = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Calculate ln(a) of the DDL layer species, according to the formula:
            // cD = c*enr,
            // where c   is the concentration of the species in the aqueous solution and
            //       enr is the enrichment factor.
            ln_a = log(aq_state.m);

            const auto& aqmix = std::any_cast<AqueousMixture>(props.extra["AqueousMixture"]);
            const auto& aqstate = std::any_cast<AqueousMixtureState>(props.extra["AqueousMixtureState"]);

            // Fetch the charges ionic species in the solution
            z_aq = aqmix.charges();

            // Fetch the stoichiometric ionic strength
            I = aqstate.Is;
        }

        if (props.extra["ComplexationSurfaceState"].has_value())
        {
            // Export surface complexation ddl_state via `extra` data member
            auto surf_state = std::any_cast<ComplexationSurfaceState>(props.extra["ComplexationSurfaceState"]);

            // Update complexation surface potential with the DDL ionic strength
            surf_state.updatePotential(I);

            // Auxiliary constant references properties of the surface
            psi = surf_state.psi;

            // Update activity coefficient using the coulombic correction factor, Appelo etal (2005), (6.47), p. 289
            ln_g = -z_aq*F*psi/(R*T); // p. 223, PHREEQC documentation (number of species in DDL = number of species in aq.sol.)
        }

        // Add the correction introduced by the activity coefficients
        ln_a += ln_g;
    };
    return fn;
}

} // namespace detail

auto ActivityModelSurfaceComplexationNoDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationNoDDL(surface_species, params);
    };
}

auto ActivityModelSurfaceComplexationWithDDL(ActivityModelSurfaceComplexationParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& surface_species)
    {
        return detail::activityModelSurfaceComplexationWithDDL(surface_species, params);
    };
}

auto ActivityModelDDL() -> ActivityModelGenerator
{
    ActivityModelDDLParams params;

    return [=](const SpeciesList& species)
    {
        return detail::activityModelDDL(species, params);
    };
}

auto ActivityModelDDL(ActivityModelDDLParams params) -> ActivityModelGenerator
{
    return [=](const SpeciesList& species)
    {
        return detail::activityModelDDL(species, params);
    };
}

} // namespace Reaktoro
