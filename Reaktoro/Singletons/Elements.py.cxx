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
#include <Reaktoro/Singletons/Elements.hpp>
using namespace Reaktoro;

void exportElements(py::module& m)
{
    py::class_<Elements, std::unique_ptr<Elements, py::nodelete>>(m, "Elements")
        .def(py::init([]() { return std::unique_ptr<Elements, py::nodelete>(&Elements::instance()); }))
        .def_static("instance", &Elements::instance, py::return_value_policy::reference)
        .def_static("data", &Elements::data, py::return_value_policy::reference)
        .def_static("clear", &Elements::clear)
        .def_static("reset", &Elements::reset)
        .def_static("append", &Elements::append)
        .def_static("replace", &Elements::replace)
        .def_static("size", &Elements::size)
        .def_static("withSymbol", &Elements::withSymbol)
        .def_static("withName", &Elements::withName)
        .def_static("withTag", &Elements::withTag)
        .def_static("withTags", &Elements::withTags)
        .def("__getitem__", [](Elements const& self, Index i) { return self.data()[i]; }, py::return_value_policy::reference)
        .def("__iter__", [](Elements const& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        ;
}
