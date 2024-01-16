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
#include <Reaktoro/Core/ElementalComposition.hpp>
using namespace Reaktoro;

void exportElementalComposition(py::module& m)
{
    py::class_<ElementalComposition>(m, "ElementalComposition")
        .def(py::init<>())
        .def(py::init<std::initializer_list<Pair<Element, double>> const&>())
        .def(py::init<Pairs<Element, double> const&>())
        .def(py::init<Pairs<String, double> const&>())
        .def("size", &ElementalComposition::size)
        .def("symbols", &ElementalComposition::symbols)
        .def("coefficients", &ElementalComposition::coefficients)
        .def("coefficient", &ElementalComposition::coefficient)
        .def("molarMass", &ElementalComposition::molarMass)
        .def("repr", &ElementalComposition::repr)
        .def("__iter__", [](ElementalComposition const& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>()) // keep object alive while iterator exists;
        ;

    py::implicitly_convertible<Pairs<Element, double>, ElementalComposition>();
    py::implicitly_convertible<Pairs<String, double>, ElementalComposition>();
}
