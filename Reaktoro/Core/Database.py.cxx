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
#include <Reaktoro/Core/Database.hpp>
using namespace Reaktoro;

void exportDatabase(py::module& m)
{
    auto addSpecies1 = [](Database& self, const Species& species)
    {
        self.addSpecies(species);
    };

    auto addSpecies2 = [](Database& self, const Vec<Species>& species)
    {
        self.addSpecies(species);
    };

    py::class_<Database>(m, "Database")
        .def(py::init<>())
        .def(py::init<const SpeciesList&>())
        .def(py::init<const ElementList&, const SpeciesList&>())
        .def(py::init<const Vec<Element>&, const Vec<Species>&>())
        .def(py::init<const Vec<Species>&>())
        .def("clear", &Database::clear)
        .def("addElement", &Database::addElement)
        .def("addSpecies", addSpecies1)
        .def("addSpecies", addSpecies2)
        .def("attachData", &Database::attachData)
        .def("extend", &Database::extend)
        .def("extendWithDatabase", &Database::extendWithDatabase)
        .def("extendWithFile", &Database::extendWithFile)
        .def("elements", &Database::elements)
        .def("species", py::overload_cast<>(&Database::species, py::const_))
        .def("speciesWithAggregateState", &Database::speciesWithAggregateState)
        .def("element", &Database::element, return_internal_ref)
        .def("species", py::overload_cast<const String&>(&Database::species, py::const_), return_internal_ref)
        .def("reaction", &Database::reaction)
        .def("attachedData", &Database::attachedData)
        .def_static("fromFile", &Database::fromFile)
        .def_static("fromEmbeddedFile", &Database::fromEmbeddedFile)
        .def_static("fromContents", &Database::fromContents)
        .def_static("fromStream", &Database::fromStream)
        .def_static("local", &Database::local)
        .def_static("embedded", &Database::embedded)
        ;
}
