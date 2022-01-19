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
#include <Reaktoro/Core/Model.py.hxx>
#include <Reaktoro/Core/ReactionThermoProps.hpp>
using namespace Reaktoro;

void exportReactionThermoProps(py::module& m)
{
    py::class_<ReactionThermoProps>(m, "ReactionThermoProps")
        .def(py::init<>())
        .def_readwrite("dG0", &ReactionThermoProps::dG0)
        .def_readwrite("dH0", &ReactionThermoProps::dH0)
        ;

    py::class_<ReactionThermoArgs>(m, "ReactionThermoArgs")
        .def_property_readonly("T", [](const ReactionThermoArgs& self) { return self.T; })
        .def_property_readonly("P", [](const ReactionThermoArgs& self) { return self.P; })
        .def_property_readonly("dV0", [](const ReactionThermoArgs& self) { return self.dV0; })
        ;

    exportModel<ReactionThermoProps, ReactionThermoArgs>(m, "ReactionThermoModel");
}
