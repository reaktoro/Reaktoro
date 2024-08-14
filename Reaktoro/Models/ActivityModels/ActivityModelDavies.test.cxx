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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelDavies.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
using namespace Reaktoro;

#define PRINT_INFO_IF_FAILS(x) INFO(#x " = \n" << std::scientific << std::setprecision(16) << x)

namespace {

/// Return mole fractions for the species.
inline auto moleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithFormula(formula); };

    ArrayXr n(species.size());
    n = 0.1;
    n[idx("H2O")] = 55.508;
    n[idx("H+" )] = 1e-7;
    n[idx("OH-")] = 1e-7;
    n[idx("Na+")] = 0.3;
    n[idx("Cl-")] = 0.3;
    n[idx("CO2")] = 0.9;

    return n / n.sum();
}

// Check if the activities of the aqueous species are correct assuming activity coefficients are.
inline auto checkActivities(ArrayXrConstRef x, ActivityPropsConstRef props)
{
    const auto iH2O = 0;

    // The concentrations of the species (molalities for solutes, mole fraction for solvent water)
    ArrayXr c = x/(x[iH2O] * waterMolarMass);
    c[iH2O] = x[iH2O];

    for(auto i = 0; i < x.size(); ++i)
    {
        INFO("i = " << i);
        CHECK( exp(props.ln_a[i] - props.ln_g[i]) == Approx(c[i]) );
    }
}

} // anonymous namespace

TEST_CASE("Testing ActivityModelDavies", "[ActivityModelDavies]")
{
    Catch::StringMaker<double>::precision = 15;

    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    SECTION("Checking the activity coefficients")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDavies()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( double(exp(props.ln_g[0]))  == Approx(0.935122413217240) ); // H2O
        CHECK( double(exp(props.ln_g[1]))  == Approx(0.800165626605902) ); // H+
        CHECK( double(exp(props.ln_g[2]))  == Approx(0.800165626605902) ); // OH-
        CHECK( double(exp(props.ln_g[3]))  == Approx(0.800165626605902) ); // Na+
        CHECK( double(exp(props.ln_g[4]))  == Approx(0.800165626605902) ); // Cl-
        CHECK( double(exp(props.ln_g[5]))  == Approx(0.409939308642971) ); // Ca++
        CHECK( double(exp(props.ln_g[6]))  == Approx(0.800165626605902) ); // HCO3-
        CHECK( double(exp(props.ln_g[7]))  == Approx(0.409939308642971) ); // CO3--
        CHECK( double(exp(props.ln_g[8]))  == Approx(1.273505728674341) ); // CO2
        CHECK( double(exp(props.ln_g[9]))  == Approx(1.273505728674341) ); // NaCl
        CHECK( double(exp(props.ln_g[10])) == Approx(1.273505728674341) ); // HCl
        CHECK( double(exp(props.ln_g[11])) == Approx(1.273505728674341) ); // NaOH

        checkActivities(x, props);
    }

    SECTION("Checking the activity coefficients using custom parameters")
    {
        // Set non-default values for the activity model
        ActivityModelDaviesParams params;
        params.bions = 0.2;
        params.bneutrals = 0.15;

        // Construct the activity props function with the given parameters and aqueous species.
        ActivityModel fn = ActivityModelDavies(params)(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( double(exp(props.ln_g[0]))  == Approx(0.923127724612841) ); // H2O
        CHECK( double(exp(props.ln_g[1]))  == Approx(0.707914308149729) ); // H+
        CHECK( double(exp(props.ln_g[2]))  == Approx(0.707914308149729) ); // OH-
        CHECK( double(exp(props.ln_g[3]))  == Approx(0.707914308149729) ); // Na+
        CHECK( double(exp(props.ln_g[4]))  == Approx(0.707914308149729) ); // Cl-
        CHECK( double(exp(props.ln_g[5]))  == Approx(0.251143973372544) ); // Ca++
        CHECK( double(exp(props.ln_g[6]))  == Approx(0.707914308149729) ); // HCO3-
        CHECK( double(exp(props.ln_g[7]))  == Approx(0.251143973372544) ); // CO3--
        CHECK( double(exp(props.ln_g[8]))  == Approx(1.437147535165123) ); // CO2
        CHECK( double(exp(props.ln_g[9]))  == Approx(1.437147535165123) ); // NaCl
        CHECK( double(exp(props.ln_g[10])) == Approx(1.437147535165123) ); // HCl
        CHECK( double(exp(props.ln_g[11])) == Approx(1.437147535165123) ); // NaOH

        checkActivities(x, props);
    }
}
