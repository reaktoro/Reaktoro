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
#include <Reaktoro/Core/Data.hpp>
using namespace Reaktoro;

void exportData(py::module& m)
{
    // Currently, Python bindings for Data is not working correctly
    // with Vec<Data>, Dict<String,Data> on clang v14.0.6 in macOS.
    // So, constructors and conversion to Data related to these types are disable
    // in the reaktoro Python package for now (in macOS only!).

    py::class_<Data>(m, "Data")
        .def(py::init<>(), "Construct a default Data instance with null value.")
        .def(py::init<bool const&>(), "Construct a Data object as a boolean value.")
        .def(py::init<String const&>(), "Construct a Data object as a string value.")
        .def(py::init<int const&>(), "Construct a Data object as a integer value.")
        .def(py::init<double const&>(), "Construct a Data object as a number value.")
        .def(py::init<real const&>(), "Construct a Data object as a number value.")
#ifndef __APPLE__
        .def(py::init<Map<String, Data> const&>(), "Construct a Data object as a dictionary object.")
        .def(py::init<Vec<Data> const&>(), "Construct a Data object as a list object.")
#endif
        .def_static("parse", py::overload_cast<Chars>(&Data::parse), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parse", py::overload_cast<String const&>(&Data::parse), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parse", py::overload_cast<std::istream&>(&Data::parse), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parseYaml", py::overload_cast<Chars>(&Data::parseYaml), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parseYaml", py::overload_cast<String const&>(&Data::parseYaml), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parseYaml", py::overload_cast<std::istream&>(&Data::parseYaml), "Return a Data object by parsing an YAML formatted string.")
        .def_static("parseJson", py::overload_cast<Chars>(&Data::parseJson), "Return a Data object by parsing a JSON formatted string.")
        .def_static("parseJson", py::overload_cast<String const&>(&Data::parseJson), "Return a Data object by parsing a JSON formatted string.")
        .def_static("parseJson", py::overload_cast<std::istream&>(&Data::parseJson), "Return a Data object by parsing a JSON formatted string.")
        .def_static("load", &Data::load, "Return a Data object by parsing either an YAML or JSON formatted file at a given path.")
        .def_static("loadYaml", &Data::loadYaml, "Return a Data object by parsing an YAML formatted file at a given path.")
        .def_static("loadJson", &Data::loadJson, "Return a Data object by parsing a JSON formatted file at a given path.")
        .def("asBoolean", &Data::asBoolean, "Return this Data object as a boolean value.")
        .def("asString", &Data::asString, return_internal_ref, "Return this Data object as a string.")
        .def("asInteger", &Data::asInteger, "Return this Data object as an integer number.")
        .def("asFloat", &Data::asFloat, "Return this Data object as a float number.")
        .def("asDict", &Data::asDict, return_internal_ref, "Return this Data object as a dictionary object.")
        .def("asList", &Data::asList, return_internal_ref, "Return this Data object as a list object.")
        .def("asNull", &Data::asNull, "Return this Data object as a nullptr value.")
        .def("isBoolean", &Data::isBoolean, "Return true if this Data object is a boolean value.")
        .def("isString", &Data::isString, "Return true if this Data object is a string.")
        .def("isInteger", &Data::isInteger, "Return true if this Data object is an integer number.")
        .def("isFloat", &Data::isFloat, "Return true if this Data object is a float number.")
        .def("isDict", &Data::isDict, "Return true if this Data object is a dictionary object.")
        .def("isList", &Data::isList, "Return true if this Data object is a list object.")
        .def("isNull", &Data::isNull, "Return true if this Data object is a null value.")
        .def("__getitem__", [](Data const& self, String const& key) { return self[key]; }, "Return the child Data object with given key, presuming this Data object is a dictionary.")
        .def("__getitem__", [](Data const& self, Index const& index) { return self[index]; }, "Return the child Data object with given index, presuming this Data object is a list.")
        .def("__getitem__", [](Data& self, String const& key) { return self[key]; }, "Return the child Data object with given key if, presuming this Data object is a dictionary.")
        .def("__getitem__", [](Data& self, Index const& index) { return self[index]; }, "Return the child Data object with given index, presuming this Data object is a list.")
        .def("at", py::overload_cast<String const&>(&Data::at, py::const_), "Return the child Data object with given key, presuming this Data object is a dictionary.")
        .def("at", py::overload_cast<Index const&>(&Data::at, py::const_), "Return the child Data object with given index, presuming this Data object is a list.")
        .def("optional", &Data::optional, "Return an optional child Data object with given key, presuming this Data object is a dictionary.")
        .def("required", &Data::required, "Return a required to exist child Data object with given key, presuming this Data object is a dictionary.")
        .def("with", &Data::with, "Return the child Data object whose `attribute` has a given `value`, presuming this Data object is a list.")
        .def("add", py::overload_cast<Data>(&Data::add), "Add a Data object to this Data object, which becomes a list if not already.")
        .def("add", py::overload_cast<String const&, Data>(&Data::add), "Add a Data object with given key to this Data object, which becomes a dictionary if not already.")
        .def("update", &Data::update, "Update this Data object with data given in another, with key-value pairs being either added or overwritten.")
        .def("reset", &Data::reset, "Reset this Data object to a null state, deleting its current stored data.")
        .def("exists", &Data::exists, "Return true if a child parameter exists with given key, presuming this Data object is a dictionary.")
        .def("dump", &Data::dump, "Return a YAML formatted string representing the state of this Data object.")
        .def("dumpYaml", &Data::dumpYaml, "Return a YAML formatted string representing the state of this Data object.")
        .def("dumpJson", &Data::dumpJson, "Return a JSON formatted string representing the state of this Data object.")
        .def("save", &Data::save, "Save the state of this Data object into a YAML formatted file.")
        .def("saveYaml", &Data::saveYaml, "Save the state of this Data object into a YAML formatted file.")
        .def("saveJson", &Data::saveJson, "Save the state of this Data object into a JSON formatted file.")
        .def("repr", &Data::repr, "Return a YAML formatted string representing the state of this Data object.")
        .def("__str__", &Data::repr, "Return a YAML formatted string representing the state of this Data object.")
        .def("__repr__", &Data::repr, "Return a YAML formatted string representing the state of this Data object.")
        ;

    py::implicitly_convertible<Chars, Data>();
    py::implicitly_convertible<String, Data>();
    py::implicitly_convertible<double, Data>();
    py::implicitly_convertible<int, Data>();
    py::implicitly_convertible<bool, Data>();
#ifndef __APPLE__
    py::implicitly_convertible<Map<String, Data>, Data>();
    py::implicitly_convertible<Vec<Data>, Data>();
#endif
}
