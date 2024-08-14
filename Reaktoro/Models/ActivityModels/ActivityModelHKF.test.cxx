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
#include <Reaktoro/Models/ActivityModels/ActivityModelHKF.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
using namespace Reaktoro;

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

TEST_CASE("Testing ActivityModelHKF", "[ActivityModelHKF]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    // Construct the activity props function with the given aqueous species.
    ActivityModel fn = ActivityModelHKF()(species);

    // Create the ActivityProps object with the results.
    ActivityProps props = ActivityProps::create(species.size());

    // Evaluate the activity props function
    fn(props, {T, P, x});

    CHECK( exp(props.ln_g[0])  == Approx(0.9979044120) ); // H2O
    CHECK( exp(props.ln_g[1])  == Approx(0.6213789582) ); // H+
    CHECK( exp(props.ln_g[2])  == Approx(0.7319869640) ); // OH-
    CHECK( exp(props.ln_g[3])  == Approx(0.6613829974) ); // Na+
    CHECK( exp(props.ln_g[4])  == Approx(0.6752174671) ); // Cl-
    CHECK( exp(props.ln_g[5])  == Approx(0.1860570160) ); // Ca++
    CHECK( exp(props.ln_g[6])  == Approx(0.6537644018) ); // HCO3-
    CHECK( exp(props.ln_g[7])  == Approx(0.1948494354) ); // CO3--
    CHECK( exp(props.ln_g[8])  == Approx(1.2735100000) ); // CO2
    CHECK( exp(props.ln_g[9])  == Approx(1.2735100000) ); // NaCl
    CHECK( exp(props.ln_g[10]) == Approx(1.2735100000) ); // HCl
    CHECK( exp(props.ln_g[11]) == Approx(1.2735100000) ); // NaOH

    checkActivities(x, props);
}
