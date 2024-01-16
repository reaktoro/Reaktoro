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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Surface.hpp>
using namespace Reaktoro;

void exportSurface(py::module& m)
{
    py::class_<Surface>(m, "Surface")
        .def(py::init<>())
        .def(py::init<String const&>())
        .def(py::init<String const&, SurfaceAreaModel const&>())
        .def("clone", &Surface::clone, "Return a deep copy of this Surface object.")
        .def("withName", &Surface::withName, "Return a duplicate of this Surface object with new name.")
        .def("withAreaModel", &Surface::withAreaModel, "Return a duplicate of this Surface object with new surface area model.")
        .def("name", &Surface::name, return_internal_ref, "Return the unique name of this surface.")
        .def("areaModel", &Surface::areaModel, return_internal_ref, "Return the area model of this surface.")
        .def("area", &Surface::area, "Calculate the area of the surface for given chemical properties of the system (in m2).")
        ;
}
