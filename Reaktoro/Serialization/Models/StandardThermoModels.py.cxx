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
#include <Reaktoro/Models/StandardThermoModels.hpp>
#include <Reaktoro/Serialization/Models/StandardThermoModels.hpp>
using namespace Reaktoro;

void exportSerializationModelsStandardThermoModels(py::module& m)
{
    //======================================================================
    // StandardThermoModelParams Types
    //======================================================================
    py::implicitly_convertible<Data, StandardThermoModelParamsConstant>();
    py::implicitly_convertible<StandardThermoModelParamsConstant, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsHKF>();
    py::implicitly_convertible<StandardThermoModelParamsHKF, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsHollandPowell>();
    py::implicitly_convertible<StandardThermoModelParamsHollandPowell, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsInterpolation>();
    py::implicitly_convertible<StandardThermoModelParamsInterpolation, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsMaierKelley>();
    py::implicitly_convertible<StandardThermoModelParamsMaierKelley, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsMineralHKF>();
    py::implicitly_convertible<StandardThermoModelParamsMineralHKF, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsNasa>();
    py::implicitly_convertible<StandardThermoModelParamsNasa, Data>();

    py::implicitly_convertible<Data, StandardThermoModelParamsWaterHKF>();
    py::implicitly_convertible<StandardThermoModelParamsWaterHKF, Data>();

    // //======================================================================
    // // ReactionStandardThermoModelParams Types
    // //======================================================================
    py::implicitly_convertible<Data, ReactionStandardThermoModelParamsConstLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsConstLgK, Data>();

    py::implicitly_convertible<Data, ReactionStandardThermoModelParamsGemsLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsGemsLgK, Data>();

    py::implicitly_convertible<Data, ReactionStandardThermoModelParamsPhreeqcLgK>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsPhreeqcLgK, Data>();

    py::implicitly_convertible<Data, ReactionStandardThermoModelParamsVantHoff>();
    py::implicitly_convertible<ReactionStandardThermoModelParamsVantHoff, Data>();

    // //======================================================================
    // // StandardVolumeModelParams Types
    // //======================================================================
    py::implicitly_convertible<Data, StandardVolumeModelParamsConstant>();
    py::implicitly_convertible<StandardVolumeModelParamsConstant, Data>();
}
