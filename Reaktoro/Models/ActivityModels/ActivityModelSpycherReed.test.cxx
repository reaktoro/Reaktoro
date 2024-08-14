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
#include <Reaktoro/Models/ActivityModels/ActivityModelSpycherReed.hpp>
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

TEST_CASE("Testing ActivityModelSpycherReed", "[ActivityModelSpycherReed]")
{
    const auto T = 300.0;
    const auto P = 12.3e5;

    WHEN("All supported gases, H2O(g), CO2(g), CH4(g), are considered.")
    {
        const auto species = SpeciesList("H2O CO2 CH4");
        const ArrayXr x = ArrayXr{{0.1, 0.6, 0.3}};

        // Construct the activity props function with the given gaseous species.
        ActivityModel fn = ActivityModelSpycherReed()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0]) == Approx(0.8867907637) ); // H2O
        CHECK( exp(props.ln_g[1]) == Approx(0.9316186544) ); // CO2
        CHECK( exp(props.ln_g[2]) == Approx(0.9966357760) ); // CH4

        checkActivities(x, P, props);
    }

    WHEN("Only H2O(g) and CO2(g) are considered.")
    {
        const auto species = SpeciesList("H2O CO2");
        const ArrayXr x = ArrayXr{{0.1, 0.9}};

        // Construct the activity props function with the given gaseous species.
        ActivityModel fn = ActivityModelSpycherReed()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0]) == Approx(0.8803341849) ); // H2O
        CHECK( exp(props.ln_g[1]) == Approx(0.9296515107) ); // CO2

        checkActivities(x, P, props);
    }

    WHEN("Only CH4(g) is considered.")
    {
        const auto species = SpeciesList("CH4");
        const ArrayXr x = ArrayXr{{1.0}};

        // Construct the activity props function with the given gaseous species.
        ActivityModel fn = ActivityModelSpycherReed()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0]) == Approx(0.9791724134) ); // CH4

        checkActivities(x, P, props);
    }
}
