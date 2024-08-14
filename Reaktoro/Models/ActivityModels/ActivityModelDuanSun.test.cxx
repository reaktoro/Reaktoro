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
#include <Reaktoro/Models/ActivityModels/ActivityModelDebyeHuckel.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelDuanSun.hpp>
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

TEST_CASE("Testing ActivityModelDuanSun", "[ActivityModelDuanSun]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    // Construct the activity props function with the given aqueous species.
    ActivityModel debyehuckelfn = ActivityModelDebyeHuckel()(species);

    // Create the ActivityProps object with the results.
    ActivityProps props = ActivityProps::create(species.size());

    // Evaluate the activity props function
    debyehuckelfn(props, {T, P, x});

    WHEN("Using ActivityModelDuanSun(CO2)")
    {
        ActivityModel fn = ActivityModelDuanSun("CO2")(species);
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.9269890137) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.7429198411) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.5772424599) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.7363279956) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.6080001197) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.2501338902) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.6538562298) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.1827801645) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.1585365027) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    WHEN("Using ActivityModelDuanSun(NaCl)")
    {
        ActivityModel fn = ActivityModelDuanSun("NaCl")(species);
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.9269890137) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.7429198411) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.5772424599) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.7363279956) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.6080001197) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.2501338902) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.6538562298) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.1827801645) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.2735057287) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.1585365027) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    WHEN("A base activity model, such as Debye-Huckel, has not been used previously")
    {
        ActivityModel fn = ActivityModelDuanSun("CO2")(species);
        props.extra = {}; // there is no AqueousMixture not AqueousMixtureState in the `extra` data member of `props`

        REQUIRE_THROWS( fn(props, {T, P, x}) );
    }
}
