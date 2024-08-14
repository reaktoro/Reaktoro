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
#include <Reaktoro/Water/WaterConstants.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ActivityModelDebyeHuckelParams", "[ActivityModelDebyeHuckel]")
{
    ActivityModelDebyeHuckelParams params;

    params.setKielland();

    REQUIRE( params.aiondefault     == 0.0 );
    REQUIRE( params.biondefault     == 0.0 );
    REQUIRE( params.bneutraldefault == 0.1 );
    REQUIRE( params.aions.size()           );
    REQUIRE( params.bions.empty()          );

    REQUIRE( params.aion("H+"  )    == 9.0 );
    REQUIRE( params.aion("Li+" )    == 6.0 );
    REQUIRE( params.aion("Rb+" )    == 2.5 );
    REQUIRE( params.aion("Cs+" )    == 2.5 );
    REQUIRE( params.aion("NH4+")    == 2.5 );
    REQUIRE( params.aion("Tl+" )    == 2.5 );
    REQUIRE( params.aion("Ag+" )    == 2.5 );
    REQUIRE( params.aion("K+"  )    == 3.0 );
    REQUIRE( params.aion("Cl-" )    == 3.0 );
    REQUIRE( params.aion("Br-" )    == 3.0 );
    REQUIRE( params.aion("I-"  )    == 3.0 );
    REQUIRE( params.aion("CN-" )    == 3.0 );

    REQUIRE( params.aion("H[+]"  )    == 9.0 );
    REQUIRE( params.aion("Li[+]" )    == 6.0 );
    REQUIRE( params.aion("Rb[+]" )    == 2.5 );
    REQUIRE( params.aion("Cs[+]" )    == 2.5 );
    REQUIRE( params.aion("NH4[+]")    == 2.5 );
    REQUIRE( params.aion("Tl[+]" )    == 2.5 );
    REQUIRE( params.aion("Ag[+]" )    == 2.5 );
    REQUIRE( params.aion("K[+]"  )    == 3.0 );
    REQUIRE( params.aion("Cl[-]" )    == 3.0 );
    REQUIRE( params.aion("Br[-]" )    == 3.0 );
    REQUIRE( params.aion("I[-]"  )    == 3.0 );
    REQUIRE( params.aion("CN[-]" )    == 3.0 );

    params.setWATEQ4F();

    REQUIRE( params.aiondefault     == 0.0 );
    REQUIRE( params.biondefault     == 0.0 );
    REQUIRE( params.bneutraldefault == 0.1 );
    REQUIRE( params.aions.size()           );
    REQUIRE( params.bions.size()           );

    REQUIRE( params.aion("Ca+2"      ) ==  5.00 );
    REQUIRE( params.aion("Mg+2"      ) ==  5.50 );
    REQUIRE( params.aion("Na+"       ) ==  4.00 );
    REQUIRE( params.aion("K+"        ) ==  3.50 );
    REQUIRE( params.aion("Cl-"       ) ==  3.50 );
    REQUIRE( params.aion("SO4-2"     ) ==  5.00 );
    REQUIRE( params.aion("HCO3-"     ) ==  5.40 );
    REQUIRE( params.aion("CO3-2"     ) ==  5.40 );
    REQUIRE( params.aion("Sr+2"      ) ==  5.26 );
    REQUIRE( params.aion("H+"        ) ==  9.00 );
    REQUIRE( params.aion("OH-"       ) ==  3.50 );
    REQUIRE( params.aion("SrHCO3+"   ) ==  5.40 );
    REQUIRE( params.aion("SrOH+"     ) ==  5.00 );
    REQUIRE( params.aion("Cu(S4)2-3" ) == 23.00 );
    REQUIRE( params.aion("CuS4S5-3"  ) == 25.00 );
    REQUIRE( params.aion("S2-2"      ) ==  6.50 );
    REQUIRE( params.aion("S3-2"      ) ==  8.00 );
    REQUIRE( params.aion("S4-2"      ) == 10.00 );
    REQUIRE( params.aion("S5-2"      ) == 12.00 );
    REQUIRE( params.aion("S6-2"      ) == 14.00 );
    REQUIRE( params.aion("Ag(S4)2-3" ) == 22.00 );
    REQUIRE( params.aion("AgS4S5-3"  ) == 24.00 );
    REQUIRE( params.aion("Ag(HS)S4-2") == 15.00 );

    REQUIRE( params.bion("Ca+2"     ) ==  0.165 );
    REQUIRE( params.bion("Mg+2"     ) ==  0.200 );
    REQUIRE( params.bion("Na+"      ) ==  0.075 );
    REQUIRE( params.bion("K+"       ) ==  0.015 );
    REQUIRE( params.bion("Cl-"      ) ==  0.015 );
    REQUIRE( params.bion("SO4-2"    ) == -0.040 );
    REQUIRE( params.bion("HCO3-"    ) ==  0.000 );
    REQUIRE( params.bion("CO3-2"    ) ==  0.000 );
    REQUIRE( params.bion("H2CO3(aq)") ==  0.000 );
    REQUIRE( params.bion("Sr+2"     ) ==  0.121 );

    params.setPHREEQC();

    REQUIRE( params.aiondefault     == 0.0 );
    REQUIRE( params.biondefault     == 0.0 );
    REQUIRE( params.bneutraldefault == 0.1 );
    REQUIRE( params.aions.size()           );
    REQUIRE( params.bions.size()           );

    REQUIRE( params.aion("Al(OH)2+" ) == 5.4  );
    REQUIRE( params.aion("Al(OH)4-" ) == 4.5  );
    REQUIRE( params.aion("Al(SO4)2-") == 4.5  );
    REQUIRE( params.aion("Al+3"     ) == 9.0  );
    REQUIRE( params.aion("AlF+2"    ) == 5.4  );
    REQUIRE( params.aion("AlF2+"    ) == 5.4  );
    REQUIRE( params.aion("AlF4-"    ) == 4.5  );
    REQUIRE( params.aion("AlOH+2"   ) == 5.4  );
    REQUIRE( params.aion("AlSO4+"   ) == 4.5  );
    REQUIRE( params.aion("Ba+2"     ) == 4.0  );
    REQUIRE( params.aion("BaOH+"    ) == 5.0  );
    REQUIRE( params.aion("Br-"      ) == 3.0  );
    REQUIRE( params.aion("CO3-2"    ) == 5.4  );
    REQUIRE( params.aion("Ca+2"     ) == 5.0  );
    REQUIRE( params.aion("CaH2PO4+" ) == 5.4  );
    REQUIRE( params.aion("CaHCO3+"  ) == 6.0  );
    REQUIRE( params.aion("CaPO4-"   ) == 5.4  );
    REQUIRE( params.aion("Cl-"      ) == 3.63 );
    REQUIRE( params.aion("Cu+"      ) == 2.5  );
    REQUIRE( params.aion("Cu+2"     ) == 6.0  );
    REQUIRE( params.aion("CuCl+"    ) == 4.0  );
    REQUIRE( params.aion("CuCl2-"   ) == 4.0  );
    REQUIRE( params.aion("CuCl3-"   ) == 4.0  );
    REQUIRE( params.aion("CuCl3-2"  ) == 5.0  );
    REQUIRE( params.aion("CuCl4-2"  ) == 5.0  );
    REQUIRE( params.aion("CuOH+"    ) == 4.0  );
    REQUIRE( params.aion("F-"       ) == 3.5  );
    REQUIRE( params.aion("Fe(OH)2+" ) == 5.4  );
    REQUIRE( params.aion("Fe(OH)3-" ) == 5.0  );
    REQUIRE( params.aion("Fe(OH)4-" ) == 5.4  );
    REQUIRE( params.aion("Fe+2"     ) == 6.0  );
    REQUIRE( params.aion("Fe+3"     ) == 9.0  );
    REQUIRE( params.aion("FeCl+2"   ) == 5.0  );
    REQUIRE( params.aion("FeCl2+"   ) == 5.0  );
    REQUIRE( params.aion("FeF+2"    ) == 5.0  );
    REQUIRE( params.aion("FeF2+"    ) == 5.0  );
    REQUIRE( params.aion("FeH2PO4+" ) == 5.4  );
    REQUIRE( params.aion("FeH2PO4+2") == 5.4  );
    REQUIRE( params.aion("FeHPO4+"  ) == 5.0  );
    REQUIRE( params.aion("FeOH+"    ) == 5.0  );
    REQUIRE( params.aion("FeOH+2"   ) == 5.0  );
    REQUIRE( params.aion("FeSO4+"   ) == 5.0  );
    REQUIRE( params.aion("H+"       ) == 9.0  );
    REQUIRE( params.aion("H2PO4-"   ) == 5.4  );
    REQUIRE( params.aion("H2SiO4-2" ) == 5.4  );
    REQUIRE( params.aion("H3SiO4-"  ) == 4.0  );
    REQUIRE( params.aion("HCO3-"    ) == 5.4  );
    REQUIRE( params.aion("HPO4-2"   ) == 5.0  );
    REQUIRE( params.aion("HS-"      ) == 3.5  );
    REQUIRE( params.aion("K+"       ) == 3.5  );
    REQUIRE( params.aion("KHPO4-"   ) == 5.4  );
    REQUIRE( params.aion("KSO4-"    ) == 5.4  );
    REQUIRE( params.aion("Li+"      ) == 6.0  );
    REQUIRE( params.aion("LiSO4-"   ) == 5.0  );
    REQUIRE( params.aion("Mg+2"     ) == 5.5  );
    REQUIRE( params.aion("MgF+"     ) == 4.5  );
    REQUIRE( params.aion("MgH2PO4+" ) == 5.4  );
    REQUIRE( params.aion("MgHCO3+"  ) == 4.0  );
    REQUIRE( params.aion("MgOH+"    ) == 6.5  );
    REQUIRE( params.aion("MgPO4-"   ) == 5.4  );
    REQUIRE( params.aion("Mn(OH)3-" ) == 5.0  );
    REQUIRE( params.aion("Mn+2"     ) == 6.0  );
    REQUIRE( params.aion("Mn+3"     ) == 9.0  );
    REQUIRE( params.aion("MnCl+"    ) == 5.0  );
    REQUIRE( params.aion("MnCl3-"   ) == 5.0  );
    REQUIRE( params.aion("MnF+"     ) == 5.0  );
    REQUIRE( params.aion("MnHCO3+"  ) == 5.0  );
    REQUIRE( params.aion("MnOH+"    ) == 5.0  );
    REQUIRE( params.aion("NH4+"     ) == 2.5  );
    REQUIRE( params.aion("NO2-"     ) == 3.0  );
    REQUIRE( params.aion("NO3-"     ) == 3.0  );
    REQUIRE( params.aion("Na+"      ) == 4.08 );
    REQUIRE( params.aion("NaHPO4-"  ) == 5.4  );
    REQUIRE( params.aion("NaSO4-"   ) == 5.4  );
    REQUIRE( params.aion("OH-"      ) == 3.5  );
    REQUIRE( params.aion("PO4-3"    ) == 4.0  );
    REQUIRE( params.aion("S-2"      ) == 5.0  );
    REQUIRE( params.aion("SO4-2"    ) == 5.0  );
    REQUIRE( params.aion("SiF6-2"   ) == 5.0  );
    REQUIRE( params.aion("Sr+2"     ) == 5.26 );
    REQUIRE( params.aion("SrHCO3+"  ) == 5.4  );
    REQUIRE( params.aion("SrOH+"    ) == 5.0  );
    REQUIRE( params.aion("Zn+2"     ) == 5.0  );
    REQUIRE( params.aion("ZnCl+"    ) == 4.0  );
    REQUIRE( params.aion("ZnCl3-"   ) == 4.0  );
    REQUIRE( params.aion("ZnCl4-2"  ) == 5.0  );

    REQUIRE( params.bion("Ba+2" ) ==  0.153 );
    REQUIRE( params.bion("Ca+2" ) ==  0.165 );
    REQUIRE( params.bion("Cl-"  ) ==  0.017 );
    REQUIRE( params.bion("K+"   ) ==  0.015 );
    REQUIRE( params.bion("Mg+2" ) ==  0.200 );
    REQUIRE( params.bion("Na+"  ) ==  0.082 );
    REQUIRE( params.bion("SO4-2") == -0.040 );
    REQUIRE( params.bion("Sr+2" ) ==  0.121 );

    REQUIRE( params.bion("Ba++" ) ==  0.153 );
    REQUIRE( params.bion("Ca++" ) ==  0.165 );
    REQUIRE( params.bion("Cl-"  ) ==  0.017 );
    REQUIRE( params.bion("K+"   ) ==  0.015 );
    REQUIRE( params.bion("Mg++" ) ==  0.200 );
    REQUIRE( params.bion("Na+"  ) ==  0.082 );
    REQUIRE( params.bion("SO4--") == -0.040 );
    REQUIRE( params.bion("Sr++" ) ==  0.121 );

    params.setLimitingLaw();

    REQUIRE( params.aiondefault     == 0.0 );
    REQUIRE( params.biondefault     == 0.0 );
    REQUIRE( params.bneutraldefault == 0.1 );
    REQUIRE( params.aions.empty()          );
    REQUIRE( params.bions.empty()          );
    REQUIRE( params.bneutrals.empty()      );
}

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

TEST_CASE("Testing ActivityModelDebyeHuckel", "[ActivityModelDebyeHuckel]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    const auto T = 300.0;
    const auto P = 12.3e5;
    const auto x = moleFractions(species);

    SECTION("Checking the activity coefficients")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDebyeHuckel()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
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
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    SECTION("Checking the activity coefficients when using PHREEQC mode")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDebyeHuckel()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
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
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    SECTION("Checking the activity coefficients when using WATEQ4F mode")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDebyeHuckelWATEQ4F()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.9265452628) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.7429198411) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.5772424599) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.7197968710) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.5985609831) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.2501338902) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.6538562298) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.1827801645) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.2735057287) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    SECTION("Checking the activity coefficients when using Kielland mode")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDebyeHuckelKielland()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.9209334298) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.7429198411) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.5772424599) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.6004257058) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.5513105763) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.2047658694) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.6004257058) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.1489680248) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.2735057287) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }

    SECTION("Checking the activity coefficients when using limiting law mode")
    {
        // Construct the activity props function with the given aqueous species.
        ActivityModel fn = ActivityModelDebyeHuckelLimitingLaw()(species);

        // Create the ActivityProps object with the results.
        ActivityProps props = ActivityProps::create(species.size());

        // Evaluate the activity props function
        fn(props, {T, P, x});

        CHECK( exp(props.ln_g[0])  == Approx(0.8469626739) ); // H2O
        CHECK( exp(props.ln_g[1])  == Approx(0.3025737114) ); // H+
        CHECK( exp(props.ln_g[2])  == Approx(0.3025737114) ); // OH-
        CHECK( exp(props.ln_g[3])  == Approx(0.3025737114) ); // Na+
        CHECK( exp(props.ln_g[4])  == Approx(0.3025737114) ); // Cl-
        CHECK( exp(props.ln_g[5])  == Approx(0.0083815583) ); // Ca++
        CHECK( exp(props.ln_g[6])  == Approx(0.3025737114) ); // HCO3-
        CHECK( exp(props.ln_g[7])  == Approx(0.0083815583) ); // CO3--
        CHECK( exp(props.ln_g[8])  == Approx(1.2735057287) ); // CO2
        CHECK( exp(props.ln_g[9])  == Approx(1.2735057287) ); // NaCl
        CHECK( exp(props.ln_g[10]) == Approx(1.2735057287) ); // HCl
        CHECK( exp(props.ln_g[11]) == Approx(1.2735057287) ); // NaOH

        checkActivities(x, props);
    }
}
