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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Kinetics/KineticsOptions.hpp>
using namespace Reaktoro;

void exportKineticsOptions(py::module& m)
{
    py::class_<KineticsOptions, EquilibriumOptions>(m, "KineticsOptions")
        .def(py::init<>())
        .def(py::init<EquilibriumOptions const&>())
        .def_readwrite("dt0", &KineticsOptions::dt0, "The time step used for preconditioning the chemical state when performing the very first chemical kinetics step.")
        ;
}
