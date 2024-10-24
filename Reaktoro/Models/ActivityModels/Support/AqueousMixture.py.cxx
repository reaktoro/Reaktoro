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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
using namespace Reaktoro;

void exportAqueousMixture(py::module& m)
{
    py::class_<AqueousMixtureState>(m, "AqueousMixtureState")
        .def(py::init<>())
        .def_readwrite("T", &AqueousMixtureState::T, "The temperature of the mixture (in K).")
        .def_readwrite("P", &AqueousMixtureState::P, "The pressure of the mixture (in Pa).")
        .def_readwrite("rho", &AqueousMixtureState::rho, "The density of water (in kg/m3)")
        .def_readwrite("epsilon", &AqueousMixtureState::epsilon, "The relative dielectric constant of water (no units)")
        .def_readwrite("Ie", &AqueousMixtureState::Ie, "The effective ionic strength of the aqueous mixture (in mol/kg)")
        .def_readwrite("Is", &AqueousMixtureState::Is, "The stoichiometric ionic strength of the aqueous mixture (in mol/kg)")
        .def_readwrite("m", &AqueousMixtureState::m, "The molalities of the aqueous species (in mol/kg)")
        .def_readwrite("ms", &AqueousMixtureState::ms, "The stoichiometric molalities of the ionic species (in mol/kg)")
        ;

    py::class_<AqueousMixture>(m, "AqueousMixture")
        .def(py::init<const SpeciesList&>())
        .def("clone", &AqueousMixture::clone, "Return a deep copy of this AqueousMixture object.")
        .def("withWaterDensityFn", &AqueousMixture::withWaterDensityFn, "Return a copy of this AqueousMixture object with replaced function for water density calculation.")
        .def("withWaterDielectricConstantFn", &AqueousMixture::withWaterDielectricConstantFn, "Return a copy of this AqueousMixture object with replaced function for water dielectric constant calculation.")
        .def("species", py::overload_cast<Index>(&AqueousMixture::species, py::const_), "Return the aqueous species in the mixture with given index.")
        .def("species", py::overload_cast<>(&AqueousMixture::species, py::const_), "Return the aqueous species in the mixture.")
        .def("neutral", &AqueousMixture::neutral, "Return the neutral aqueous solutes in the mixture.")
        .def("charged", &AqueousMixture::charged, "Return the charged aqueous solutes in the mixture.")
        .def("cations", &AqueousMixture::cations, "Return the cation solutes in the mixture.")
        .def("anions", &AqueousMixture::anions, "Return the anion solutes in the mixture.")
        .def("water", &AqueousMixture::water, "Return the aqueous solvent species in the mixture.")
        .def("indicesNeutral", &AqueousMixture::indicesNeutral, "Return the indices of the neutral aqueous solutes in the mixture.")
        .def("indicesCharged", &AqueousMixture::indicesCharged, "Return the indices of the charged aqueous solutes in the mixture.")
        .def("indicesCations", &AqueousMixture::indicesCations, "Return the indices of the cations in the mixture.")
        .def("indicesAnions", &AqueousMixture::indicesAnions, "Return the indices of the anions in the mixture.")
        .def("indexWater", &AqueousMixture::indexWater, "Return the index of the solvent species in the mixture.")
        .def("charges", &AqueousMixture::charges, "Return the electric charges of the aqueous species in the mixture.")
        .def("dissociationMatrix", &AqueousMixture::dissociationMatrix, "Return the dissociation matrix of the neutral species into charged species.")
        .def("state", [](AqueousMixture& self, real T, real P, ArrayXrConstRef x) { return self.state(T, P, x); }, "Calculate the state of the aqueous mixture.")
        .def("state", [](AqueousMixture& self, real T, real P, py::array_t<double> const& x) { return self.state(T, P, ArrayXr(ArrayXd::Map(x.data(), x.size()))); }, "Calculate the state of the aqueous mixture.")
        .def("state", &AqueousMixture::state, "Calculate the state of the aqueous mixture.")
        .def_static("setDefaultWaterDensityFn", AqueousMixture::setDefaultWaterDensityFn, "Set the default function for water density calculation when creating AqueousMixture objects.")
        .def_static("setDefaultWaterDielectricConstantFn", AqueousMixture::setDefaultWaterDielectricConstantFn, "Set the default function for water dielectric constant calculation when creating AqueousMixture objects.")
        .def_static("resetDefaultWaterDensityFn", [](){ AqueousMixture::resetDefaultWaterDensityFn(); }, "Reset the default function for water density calculation when creating AqueousMixture objects.")
        .def_static("resetDefaultWaterDielectricConstantFn", [](){ AqueousMixture::resetDefaultWaterDielectricConstantFn(); }, "Reset the default function for water dielectric constant calculation when creating AqueousMixture objects.")
        ;
}
