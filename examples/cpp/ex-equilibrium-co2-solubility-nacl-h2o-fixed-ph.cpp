// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021
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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    // Initialize a thermodynamic database
    Database db("supcrtbl.yaml");

    // Define aqueous phase with selected species
    AqueousPhase aqueousphase("H2O(aq) CO2(aq) CO3-2 Cl- H+ H2(aq) HCO3- Na+ NaCl(aq) NaOH(aq) O2(aq) OH- HCl(aq)");
    aqueousphase.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond("CO2")
    ));

    // Create a gaseous phase
    GaseousPhase gaseousphase("CO2(g) H2O(g)");
    gaseousphase.setActivityModel(ActivityModelPengRobinson());

    // Collecting all above-defined phases
    Phases phases(db);
    phases.add(aqueousphase);
    phases.add(gaseousphase);

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Specify conditions to be satisfied at chemical equilibrium
    EquilibriumSpecs specs(system);
    specs.temperature();
    specs.pressure();
    specs.pH();

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(specs);

    // Define conditions to be satisfied at chemical equilibrium
    EquilibriumConditions conditions(specs);
    conditions.temperature(60.0, "celsius");
    conditions.pressure(100.0, "bar");
    conditions.pH(4.0);
    // Define initial species amounts
    conditions.startWith("H2O(aq)", 1.0, "kg");
    conditions.startWith("Na+",     1.0, "mol");
    conditions.startWith("Cl-",     1.0, "mol");
    conditions.startWith("CO2(g)", 10.0, "mol");

    // Define initial equilibrium state
    ChemicalState state(system);

    // Equilibrate given initial state
    solver.solve(state, conditions);

    // Obtain species composition from the equilibrated state
    const auto n = state.speciesAmounts();

    // Print the species and theirs amounts
    std::cout << std::setw(20) << "Species"
              << std::setw(20) << "Amount" << std::endl;

    for(auto i = 0; i < n.size(); ++i)
    {
        // Print only species with nonzero amounts
        if(n[i] > 1e-16)
            std::cout << std::setw(20) << system.species(i).name()
                      << std::setw(20) << n[i] << std::endl;
    }

    std::cout << state << std::endl;

    return 0;
}