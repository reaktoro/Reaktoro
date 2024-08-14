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
#include <Reaktoro/catch.hxx>

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelCubicEOS.hpp>
using namespace Reaktoro;

namespace {

// Check if the activities of the fluid species are correct assuming activity coefficients are.
inline auto checkActivities(ArrayXrConstRef x, real P, ActivityPropsConstRef props)
{
    const auto Pbar = P * 1e-5;

    // The concentrations of the species (partial pressures in bar)
    ArrayXr c = x * Pbar;

    for(auto i = 0; i < x.size(); ++i)
    {
        INFO("i = " << i);
        CHECK( exp(props.ln_a[i] - props.ln_g[i]) == Approx(c[i]) );
    }
}

} // anonymous namespace

TEST_CASE("Testing ActivityModelCubicEOS", "[ActivityModelCubicEOS]")
{
    const auto T = 300.0;
    const auto P = 12.3e5;

    //=============================================
    // FLUID MIXTURE (CO2, H2O, CH4)
    //=============================================
    WHEN("The gases are CO2, H2O, CH4")
    {
        const auto species = SpeciesList("CO2 H2O CH4");

        const ArrayXr x = ArrayXr{{0.90, 0.08, 0.02}};

        // Construct the activity props function with the given species.
        ActivityModel fn = ActivityModelPengRobinson()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        WHEN("Conditions correspond to gas state")
        {
            const auto T = 25.0 + 273.15; // 25 °C
            const auto P = 1.0 * 1e5;     // 1 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(0.02464) );

            CHECK( exp(props.ln_g[0]) == Approx(0.994555) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.986592) ); // H2O
            CHECK( exp(props.ln_g[2]) == Approx(0.998501) ); // CH4

            CHECK( props.som == StateOfMatter::Gas );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to liquid state")
        {
            const auto T = 10.0 + 273.15; // 10 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(4.21022e-05) );

            CHECK( exp(props.ln_g[0]) == Approx(0.375122) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.00662715) ); // H2O
            CHECK( exp(props.ln_g[2]) == Approx(2.21312) ); // CH4

            CHECK( props.som == StateOfMatter::Liquid );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to supercritical state")
        {
            const auto T = 60.0 + 273.15; // 60 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(9.33293e-05) );

            CHECK( exp(props.ln_g[0]) == Approx(0.676616) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.137854) ); // H2O
            CHECK( exp(props.ln_g[2]) == Approx(1.34433) ); // CH4

            CHECK( props.som == StateOfMatter::Supercritical );

            checkActivities(x, P, props);
        }
    }

    //=============================================
    // FLUID MIXTURE (CO2, H2O)
    //=============================================
    WHEN("The gases are CO2 and H2O")
    {
        const auto species = SpeciesList("CO2 H2O");

        const ArrayXr x = ArrayXr{{0.9, 0.1}};

        // Construct the activity props function with the given species.
        ActivityModel fn = ActivityModelPengRobinson()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        WHEN("Conditions correspond to gas state")
        {
            const auto T = 25.0 + 273.15; // 25 °C
            const auto P = 1.0 * 1e5;     // 1 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(0.024634) );

            CHECK( exp(props.ln_g[0]) == Approx(0.994566) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.986475) ); // H2O

            CHECK( props.som == StateOfMatter::Gas );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to liquid state")
        {
            const auto T = 10.0 + 273.15; // 10 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(4.04694e-05) );

            CHECK( exp(props.ln_g[0]) == Approx(0.382635) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.00536954) ); // H2O

            CHECK( props.som == StateOfMatter::Liquid );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to supercritical state")
        {
            const auto T = 60.0 + 273.15; // 60 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(6.92881e-05) );

            CHECK( exp(props.ln_g[0]) == Approx(0.707373) ); // CO2
            CHECK( exp(props.ln_g[1]) == Approx(0.0855991) ); // H2O

            CHECK( props.som == StateOfMatter::Supercritical );

            checkActivities(x, P, props);
        }
    }

    //=============================================
    // FLUID MIXTURE (CO2)
    //=============================================
    WHEN("There is only CO2")
    {
        const auto species = SpeciesList("CO2");

        const ArrayXr x = ArrayXr{{1.0}};

        // Construct the activity props function with the given species.
        ActivityModel fn = ActivityModelPengRobinson()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        WHEN("Conditions correspond to gas state")
        {
            const auto T = 25.0 + 273.15; // 25 °C
            const auto P = 1.0 * 1e5;     // 1 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(0.0246538) );

            CHECK( exp(props.ln_g[0]) == Approx(0.994544) ); // CO2

            CHECK( props.som == StateOfMatter::Gas );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to liquid state")
        {
            const auto T = 10.0 + 273.15; // 10 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(4.80345e-05) );

            CHECK( exp(props.ln_g[0]) == Approx(0.36694) ); // CO2

            CHECK( props.som == StateOfMatter::Liquid );

            checkActivities(x, P, props);
        }

        WHEN("Conditions correspond to supercritical state")
        {
            const auto T = 60.0 + 273.15; // 60 °C
            const auto P = 100.0 * 1e5;   // 100 bar

            fn(props, {T, P, x}); // Evaluate the activity props function

            CHECK( props.Vx == Approx(0.000150429) );

            CHECK( exp(props.ln_g[0]) == Approx(0.657945) ); // CO2

            CHECK( props.som == StateOfMatter::Supercritical );

            checkActivities(x, P, props);
        }
    }
}
