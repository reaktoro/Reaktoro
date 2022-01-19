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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
using namespace Reaktoro;

void exportEquilibriumProjector(py::module& m)
{
    py::class_<EquilibriumPredictor>(m, "EquilibriumPredictor")
        .def(py::init<const ChemicalState&, const EquilibriumSensitivity&>())
        .def("predict", py::overload_cast<ChemicalState&, const EquilibriumConditions&>(&EquilibriumPredictor::predict))
        .def("predict", py::overload_cast<ChemicalState&, const EquilibriumConditions&, VectorXrConstRef>(&EquilibriumPredictor::predict))
        ;
}
