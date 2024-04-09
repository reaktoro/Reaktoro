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

#pragma once

//----------------------------------------------------------------------------------------------------------------
// Add defines above to allow pybind11 packages produced with different compilers/versions to be used together.
// https://github.com/pybind/pybind11/pull/2602
//----------------------------------------------------------------------------------------------------------------
#define PYBIND11_COMPILER_TYPE ""
#define PYBIND11_STDLIB ""
#define PYBIND11_BUILD_ABI ""
//----------------------------------------------------------------------------------------------------------------

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/cast.h>
namespace py = pybind11;
using namespace py::literals;

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

const auto return_internal_ref = py::return_value_policy::reference_internal;

PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXr);
PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXrRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::ArrayXrConstRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXr);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXrRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::VectorXrConstRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::MatrixXr);
PYBIND11_MAKE_OPAQUE(Reaktoro::MatrixXrRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::MatrixXrConstRef);
PYBIND11_MAKE_OPAQUE(Reaktoro::Indices);
