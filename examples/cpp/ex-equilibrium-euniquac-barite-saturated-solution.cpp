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

// -----------------------------------------------------------------------------
// 👏 Acknowledgements 👏
// -----------------------------------------------------------------------------
// This example was originally authored by:
//   • Allan Leal (19 October 2023)
// -----------------------------------------------------------------------------

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    const auto db = Database::embedded("ExtendedUNIQUAC.v2024.yaml");
    const auto params = Params::embedded("ExtendedUNIQUAC.v2024.yaml");

    AqueousPhase solution(speciate("H O Na Ba Cl C S"));
    solution.setActivityModel(ActivityModelExtendedUNIQUAC(params));

    ChemicalSystem system(db, solution);

    ChemicalState state(system);
    state.temperature(60.0, "celsius");
    state.pressure(100.0, "bar");
    state.set("H2O"  , 1.0, "kg");
    state.set("Na+"  , 4.0, "mol");
    state.set("Cl-"  , 3.0, "mol");
    state.set("Ba+2" , 0.5, "mol");
    state.set("SO4-2", 1.0, "mol");
    state.set("CO2"  , 0.2, "mol");

    const auto result = equilibrate(state);

    AqueousProps aprops(state);
    aprops.output("aprops.txt");

    std::cout << "Check the generated file `aprops.txt` for the properties of the aqueous solution." << std::endl;
}
