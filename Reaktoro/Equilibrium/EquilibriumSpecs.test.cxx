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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSpecs", "[EquilibriumSpecs]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumSpecs specs(system);

    ChemicalState state(system);
    state.setTemperature(50.0, "celsius");
    state.setPressure(100.0, "bar");
    state.setSpeciesAmounts(1.0);

    ChemicalProps props(state);

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        CHECK( specs.numInputs()                             == 2 ); // T, P
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 0 );
        CHECK( specs.numControlVariablesP()                  == 0 );
        CHECK( specs.numControlVariablesQ()                  == 0 );
        CHECK( specs.numTitrants()                           == 0 );
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 0 );
        CHECK( specs.numConstraints()                        == 0 );
        CHECK( specs.namesInputs()                           == Strings{"T", "P"} );
        CHECK( specs.isTemperatureUnknown()                  == false );
        CHECK( specs.isPressureUnknown()                     == false );
        CHECK( specs.indexControlVariableTemperature()       == Index(-1) );
        CHECK( specs.indexControlVariablePressure()          == Index(-1) );
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        CHECK( specs.numInputs()                             == 2 ); // T, V
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 1 ); // P
        CHECK( specs.numControlVariablesP()                  == 1 ); // P
        CHECK( specs.numControlVariablesQ()                  == 0 );
        CHECK( specs.numTitrants()                           == 0 );
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 0 );
        CHECK( specs.numConstraints()                        == 1 ); // V = V(given)
        CHECK( specs.namesInputs()                           == Strings{"T", "V"} );
        CHECK( specs.namesControlVariables()                 == Strings{"P"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume"} );
        CHECK( specs.isTemperatureUnknown()                  == false );
        CHECK( specs.isPressureUnknown()                     == true );
        CHECK( specs.indexControlVariableTemperature()       == Index(-1) );
        CHECK( specs.indexControlVariablePressure()          == 0 );
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        CHECK( specs.numInputs()                             == 2 ); // V, U
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 2 ); // T, P
        CHECK( specs.numControlVariablesP()                  == 2 ); // T, P
        CHECK( specs.numControlVariablesQ()                  == 0 );
        CHECK( specs.numTitrants()                           == 0 );
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 0 );
        CHECK( specs.numConstraints()                        == 2 ); // V = V(given), U = U(given)
        CHECK( specs.namesInputs()                           == Strings{"V", "U"} );
        CHECK( specs.namesControlVariables()                 == Strings{"T", "P"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume", "internalEnergy"} );
        CHECK( specs.isTemperatureUnknown()                  == true );
        CHECK( specs.isPressureUnknown()                     == true );
        CHECK( specs.indexControlVariableTemperature()       == 0 );
        CHECK( specs.indexControlVariablePressure()          == 1 );
    }

    WHEN("temperature, pressure, and pH are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        CHECK( specs.numInputs()                             == 3 ); // T, P, pH
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 1 ); // n[H+]
        CHECK( specs.numControlVariablesP()                  == 0 );
        CHECK( specs.numControlVariablesQ()                  == 1 ); // n[H+]
        CHECK( specs.numTitrants()                           == 1 ); // [H+]
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 1 ); // [H+]
        CHECK( specs.numConstraints()                        == 1 ); // pH = pH(given)
        CHECK( specs.namesInputs()                           == Strings{"T", "P", "pH"} );
        CHECK( specs.namesControlVariables()                 == Strings{"[H+]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[H+]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[H+]"} );
        CHECK( specs.namesConstraints()                      == Strings{"pH"} );
        CHECK( specs.isTemperatureUnknown()                  == false );
        CHECK( specs.isPressureUnknown()                     == false );
        CHECK( specs.indexControlVariableTemperature()       == Index(-1) );
        CHECK( specs.indexControlVariablePressure()          == Index(-1) );
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input variables")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        CHECK( specs.numInputs()                             == 3 ); // V, S, a[CO2(g)]
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 3 ); // T, P, n[CO2]
        CHECK( specs.numControlVariablesP()                  == 2 ); // T, P
        CHECK( specs.numControlVariablesQ()                  == 1 ); // n[CO2]
        CHECK( specs.numTitrants()                           == 1 ); // [CO2]
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 1 ); // [CO2]
        CHECK( specs.numConstraints()                        == 3 ); // V = V(given), S = S(given), a[CO2(g)] = a[CO2(g)](given)
        CHECK( specs.namesInputs()                           == Strings{"V", "S", "ln(a[CO2(g)])"} );
        CHECK( specs.namesControlVariables()                 == Strings{"T", "P", "[CO2(g)]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[CO2]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[CO2]"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume", "entropy", "ln(a[CO2(g)])"} );
        CHECK( specs.isTemperatureUnknown()                  == true );
        CHECK( specs.isPressureUnknown()                     == true );
        CHECK( specs.indexControlVariableTemperature()       == 0 );
        CHECK( specs.indexControlVariablePressure()          == 1 );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");

        CHECK( specs.numInputs()                             == 6 ); // T, P, V, U, pH, pE
        CHECK( specs.numParams()                             == 0 );
        CHECK( specs.numControlVariables()                   == 4 ); // n[CO2], n[CH4], n[H+], n[e-]
        CHECK( specs.numControlVariablesP()                  == 2 ); // n[CO2], n[CH4]
        CHECK( specs.numControlVariablesQ()                  == 2 ); // n[H+], n[e-]
        CHECK( specs.numTitrants()                           == 4 ); // [CO2], [CH4], [H+], [e-]
        CHECK( specs.numTitrantsExplicit()                   == 2 ); // [CO2], [CH4]
        CHECK( specs.numTitrantsImplicit()                   == 2 ); // [H+], [e-]
        CHECK( specs.numConstraints()                        == 4 ); // V = V(given), U = U(given), pH = pH(given), pE = pE(given)
        CHECK( specs.namesInputs()                           == Strings{"T", "P", "V", "U", "pH", "pE"} );
        CHECK( specs.namesControlVariables()                 == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrantsExplicit()                 == Strings{"[CO2]", "[CH4]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[H+]", "[e-]"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume", "internalEnergy", "pH", "pE"} );
        CHECK( specs.isTemperatureUnknown()                  == false );
        CHECK( specs.isPressureUnknown()                     == false );
        CHECK( specs.indexControlVariableTemperature()       == Index(-1) );
        CHECK( specs.indexControlVariablePressure()          == Index(-1) );
    }

    WHEN("model parameters are among the input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.addInput("V");
        specs.addInput(Param("G0[H2O]", 1.0));

        CHECK( specs.numInputs()                             == 4 ); // T, P, V, G0[H2O]
        CHECK( specs.numParams()                             == 1 ); // G0[H2O]
        CHECK( specs.numControlVariables()                   == 0 );
        CHECK( specs.numControlVariablesP()                  == 0 );
        CHECK( specs.numControlVariablesQ()                  == 0 );
        CHECK( specs.numTitrants()                           == 0 );
        CHECK( specs.numTitrantsExplicit()                   == 0 );
        CHECK( specs.numTitrantsImplicit()                   == 0 );
        CHECK( specs.numConstraints()                        == 0 );
        CHECK( specs.namesInputs()                           == Strings{"T", "P", "V", "G0[H2O]"} );
        CHECK( specs.namesParams()                           == Strings{"G0[H2O]"} );
        CHECK( specs.namesControlVariables()                 == Strings{} );
        CHECK( specs.namesTitrants()                         == Strings{} );
        CHECK( specs.namesTitrantsExplicit()                 == Strings{} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{} );
        CHECK( specs.namesConstraints()                      == Strings{} );
        CHECK( specs.indicesParams()                         == Indices{3} ); // index of G0[H2O]
        CHECK( specs.isTemperatureUnknown()                  == false );
        CHECK( specs.isPressureUnknown()                     == false );
        CHECK( specs.indexControlVariableTemperature()       == Index(-1) );
        CHECK( specs.indexControlVariablePressure()          == Index(-1) );
    }

    SECTION("Checking lambda functions in equation constraints")
    {
        const auto& econstraints = specs.constraintsEquationType();

        specs.volume();
        specs.internalEnergy();
        specs.enthalpy();
        specs.gibbsEnergy();
        specs.helmholtzEnergy();
        specs.entropy();

        const VectorXr w = random(6);

        CHECK( econstraints[0].fn(props, w) == props.volume() - w[0] );
        CHECK( econstraints[1].fn(props, w) == props.internalEnergy() - w[1] );
        CHECK( econstraints[2].fn(props, w) == props.enthalpy() - w[2] );
        CHECK( econstraints[3].fn(props, w) == props.gibbsEnergy() - w[3] );
        CHECK( econstraints[4].fn(props, w) == props.helmholtzEnergy() - w[4] );
        CHECK( econstraints[5].fn(props, w) == props.entropy() - w[5] );
    }

    SECTION("Checking lambda functions in chemical potential constraints")
    {
        const auto constraintEh = GENERATE(true, false); // constrain Eh if true, pE if false (cannot be both constrained because of same titrant e-)

        specs.chemicalPotential("H2O(aq)");
        specs.lnActivity("CH4(g)");
        specs.lgActivity("CO2(g)");
        specs.activity("Ca++(aq)");
        specs.fugacity("O2");
        specs.pH();
        specs.pMg();

        if(constraintEh) specs.Eh(); else specs.pE(); // Not possible to constraint pE and Eh simultaneously

        const VectorXr w = random(8).cwiseAbs(); // the input variables for the constrained properties above

        const auto T  = props.temperature();
        const auto P  = props.pressure();
        const auto RT = universalGasConstant * T;
        const auto F  = faradayConstant;

        const auto u0CH4  = system.species().get("CH4(g)").standardThermoProps(T, P).G0;
        const auto u0CO2  = system.species().get("CO2(g)").standardThermoProps(T, P).G0;
        const auto u0Capp = system.species().get("Ca++(aq)").standardThermoProps(T, P).G0;
        const auto u0O2   = system.species().get("O2(g)").standardThermoProps(T, P).G0;
        const auto u0Hp   = system.species().get("H+(aq)").standardThermoProps(T, P).G0;
        const auto u0Mgpp = system.species().get("Mg++(aq)").standardThermoProps(T, P).G0;

        const auto& qvars = specs.controlVariablesQ();

        CHECK( qvars[0].fn(props, w) == Approx(w[0]) );
        CHECK( qvars[1].fn(props, w) == Approx(u0CH4 + RT*w[1]) );
        CHECK( qvars[2].fn(props, w) == Approx(u0CO2 + RT*w[2]) );
        CHECK( qvars[3].fn(props, w) == Approx(u0Capp + RT*w[3]) );
        CHECK( qvars[4].fn(props, w) == Approx(u0O2 + RT*log(w[4])) );
        CHECK( qvars[5].fn(props, w) == Approx(u0Hp + RT*w[5] * (-ln10)) );
        CHECK( qvars[6].fn(props, w) == Approx(u0Mgpp + RT*w[6] * (-ln10)) );
        CHECK( qvars[7].fn(props, w) == Approx(constraintEh ? -F * w[7] : RT*w[7] * (-ln10)) );
    }

    SECTION("Checking when chemical potential unknowns are introduced")
    {
        WHEN("temperature and pressure are given inputs")
        {
            specs.temperature();
            specs.pressure();

            specs.addUnknownActivity("Ca++(aq)");
            specs.addUnknownActivityCoefficient("H+(aq)");
            specs.addUnknownChemicalPotential("CaCO3(s)");
            specs.addUnknownStandardChemicalPotential("SiO2(s)");

            CHECK( specs.numControlVariablesP() == 4 ); // one for each addUnknown* call above

            auto pvars = specs.controlVariablesP();

            CHECK( pvars[0].ispecies == system.species().index("Ca++(aq)") );
            CHECK( pvars[1].ispecies == system.species().index("H+(aq)") );
            CHECK( pvars[2].ispecies == system.species().index("CaCO3(s)") );
            CHECK( pvars[3].ispecies == system.species().index("SiO2(s)") );

            CHECK( pvars[0].name == "a[Ca++(aq)]" );
            CHECK( pvars[1].name == "g[H+(aq)]" );
            CHECK( pvars[2].name == "u[CaCO3(s)]" );
            CHECK( pvars[3].name == "u0[SiO2(s)]" );

            const auto pk = 3.7;
            const auto ln_pk = log(pk);

            const auto T = props.temperature();
            const auto P = props.pressure();
            const auto RT = universalGasConstant * T;

            const auto species0 = system.species(pvars[0].ispecies);
            const auto species1 = system.species(pvars[1].ispecies);
            const auto species2 = system.species(pvars[2].ispecies);
            const auto species3 = system.species(pvars[3].ispecies);

            const auto G0_0 = species0.standardThermoProps(T, P).G0;
            const auto G0_1 = species1.standardThermoProps(T, P).G0;
            const auto G0_2 = species2.standardThermoProps(T, P).G0;
            const auto G0_3 = species3.standardThermoProps(T, P).G0;

            const auto ln_c1 = props.speciesConcentrationLn(pvars[1].ispecies);
            const auto ln_a3 = props.speciesActivityLn(pvars[3].ispecies);

            CHECK( pvars[0].fn(props, pk) == Approx(G0_0 + RT*log(pk)) );
            CHECK( pvars[1].fn(props, pk) == Approx(G0_1 + RT*(log(pk) + ln_c1)) );
            CHECK( pvars[2].fn(props, pk) == Approx(pk) );
            CHECK( pvars[3].fn(props, pk) == Approx(pk + RT*ln_a3) );
        }
    }
}
