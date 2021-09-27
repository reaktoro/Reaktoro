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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Svetlana Kyas (27 Spetember 2021)
//
// and since revised by:
//   •
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Initialize a thermodynamic database
    Database db("supcrtbl.yaml");

    // Create an aqueous phase automatically selecting all species with given elements, excluding species with tag `organic`
    AqueousPhase aqueousphase("Ca(HCO3)+ CO3-2 CO2(aq) CaCO3(aq) Ca+2 CaSO4(aq) CaOH+ Cl- "
                              "FeCl+2 FeCl2(aq) FeCl+ Fe+2 FeOH+ FeOH+2 Fe+3 "
                              "H2(aq) HSO4- H2S(aq) HS- H2O(aq) H+ HCO3- K+ KSO4- "
                              "Mg+2 MgSO4(aq) MgCO3(aq) MgOH+ Mg(HCO3)+ Na+ NaSO4- "
                              "O2(aq) OH- S5-2 S4-2 S3-2 S2-2 SO4-2");
    aqueousphase.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond("CO2")
    ));

    // Create a mineral phases
    MineralPhases mineralphases("Siderite Pyrite Hematite");

    // Collect all above-defined phases
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(mineralphases);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();
    specs.pE();

    // Define equilibrium solver
    EquilibriumSolver solver(specs);

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(25.0, "celsius");
    conditions.pressure(1.0, "atm");
    conditions.pH(5.726);
    conditions.pE(8.220);

    // Specify chemical composition to start from
    conditions.startWith("H2O(aq)" , 58.0     , "kg");
    conditions.startWith("Cl-"     , 1122.3e-3, "kg");
    conditions.startWith("Na+"     , 624.08e-3, "kg");
    conditions.startWith("SO4-2"   , 157.18e-3, "kg");
    conditions.startWith("Mg+2"    , 74.820e-3, "kg");
    conditions.startWith("Ca+2"    , 23.838e-3, "kg");
    conditions.startWith("K+"      , 23.142e-3, "kg");
    conditions.startWith("HCO3-"   , 8.236e-3 , "kg");
    conditions.startWith("O2(aq)"  , 58e-12   , "kg");
    conditions.startWith("Siderite", 0.0      , "mol");
    conditions.startWith("Pyrite"  , 0.0      , "mol");
    conditions.startWith("Hematite", 0.0      , "mol");
    conditions.startWith("HS-"     , 0.0196504, "mol");
    conditions.startWith("H2S(aq)" , 0.167794 , "mol");

    // Create chemical state
    ChemicalState state(system);

    // Equilibrate the initial state with equilibrium condition
    EquilibriumResult result = solver.solve(state, conditions);

    // Output the chemical state to a text file
    state.output("state.txt");

    // Output characteristics of equilibrium calculations
    std::cout << "Equilibrium calculation characteristics: " << std::endl;
    std::cout << " * iterations = " << result.optima.iterations << std::endl;
    std::cout << " * succeeded  = " << result.optima.succeeded << std::endl;

    // Compute aqueous properties at the calculated equilibrium state
    AqueousProps aprops(state);

    // Output the pH, pE, and ionic strength of the solution after equilibration
    std::cout << "pH after equilibration:             " << aprops.pH() << std::endl;
    std::cout << "pE after equilibration:             " << aprops.pE() << std::endl;
    std::cout << "Ionic strength after equilibration: " << aprops.ionicStrength() << std::endl;

    return 0;
}