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
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
using namespace Reaktoro;

template<typename ChemicalPropsPhaseType>
auto createTemplateClassForChemicalPropsPhaseType(py::module& m, const char* name)
{
    return py::class_<ChemicalPropsPhaseType>(m, name)
        .def("phase", &ChemicalPropsPhaseType::phase, return_internal_ref, "Return the underlying Phase object.")
        .def("data", &ChemicalPropsPhaseType::data, return_internal_ref, "Return the primary chemical property data of the phase from which others are calculated.")
        .def("stateOfMatter", &ChemicalPropsPhaseType::stateOfMatter, "Return the state of matter of the phase.")
        .def("temperature", &ChemicalPropsPhaseType::temperature, "Return the temperature of the phase (in K).")
        .def("pressure", &ChemicalPropsPhaseType::pressure, "Return the pressure of the phase (in Pa).")
        .def("speciesAmounts", &ChemicalPropsPhaseType::speciesAmounts, return_internal_ref, "Return the amounts of the species in the phase (in mol).")
        .def("speciesMoleFractions", &ChemicalPropsPhaseType::speciesMoleFractions, return_internal_ref, "Return the mole fractions of the species in the phase.")
        .def("speciesActivityCoefficientsLn", &ChemicalPropsPhaseType::speciesActivityCoefficientsLn, return_internal_ref, "Return the ln activity coefficients of the species in the phase.")
        .def("speciesActivitiesLn", &ChemicalPropsPhaseType::speciesActivitiesLn, return_internal_ref, "Return the ln activities of the species in the phase.")
        .def("speciesChemicalPotentials", &ChemicalPropsPhaseType::speciesChemicalPotentials, return_internal_ref, "Return the chemical potentials of the species in the phase (in J/mol).")
        .def("speciesPartialMolarVolumes", &ChemicalPropsPhaseType::speciesPartialMolarVolumes, "Return the partial molar volumes of the species in the phase (in m³/mol).")
        .def("speciesStandardVolumes", &ChemicalPropsPhaseType::speciesStandardVolumes, return_internal_ref, "Return the standard partial molar volumes of the species in the phase (in m³/mol).")
        .def("speciesStandardVolumesT", &ChemicalPropsPhaseType::speciesStandardVolumesT, return_internal_ref, "Return the temperature derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).")
        .def("speciesStandardVolumesP", &ChemicalPropsPhaseType::speciesStandardVolumesP, return_internal_ref, "Return the pressure derivative of the standard molar volumes of the species in the phase (in m³/(mol·K)).")
        .def("speciesStandardGibbsEnergies", &ChemicalPropsPhaseType::speciesStandardGibbsEnergies, return_internal_ref, "Return the standard partial molar Gibbs energies of formation of the species in the phase (in J/mol).")
        .def("speciesStandardEnthalpies", &ChemicalPropsPhaseType::speciesStandardEnthalpies, return_internal_ref, "Return the standard partial molar enthalpies of formation of the species in the phase (in J/mol).")
        .def("speciesStandardEntropies", &ChemicalPropsPhaseType::speciesStandardEntropies, "Return the standard partial molar entropies of formation of the species in the phase (in J/(mol·K)).")
        .def("speciesStandardInternalEnergies", &ChemicalPropsPhaseType::speciesStandardInternalEnergies, "Return the standard partial molar internal energies of formation of the species in the phase (in J/mol).")
        .def("speciesStandardHelmholtzEnergies", &ChemicalPropsPhaseType::speciesStandardHelmholtzEnergies, "Return the standard partial molar Helmholtz energies of formation of the species in the phase (in J/mol).")
        .def("speciesStandardHeatCapacitiesConstP", &ChemicalPropsPhaseType::speciesStandardHeatCapacitiesConstP, return_internal_ref, "Return the standard partial molar isobaric heat capacities of the species in the phase (in J/(mol·K)).")
        .def("speciesStandardHeatCapacitiesConstV", &ChemicalPropsPhaseType::speciesStandardHeatCapacitiesConstV, return_internal_ref, "Return the standard partial molar isochoric heat capacities of the species in the phase (in J/(mol·K)).")
        .def("molarMass", &ChemicalPropsPhaseType::molarMass, "Return the molar mass of the phase (in kg/mol).")
        .def("molarVolume", &ChemicalPropsPhaseType::molarVolume, "Return the molar volume of the phase (in m³/mol).")
        .def("molarVolumeT", &ChemicalPropsPhaseType::molarVolumeT, "Return the temperature derivative of the molar volume of the phase (in m³/(mol·K)).")
        .def("molarVolumeP", &ChemicalPropsPhaseType::molarVolumeP, "Return the pressure derivative of the molar volume of the phase (in m³/(mol·Pa)).")
        .def("correctiveMolarVolume", &ChemicalPropsPhaseType::correctiveMolarVolume, "Return the corrective molar volume of the phase (in m³/mol).")
        .def("correctiveMolarVolumeT", &ChemicalPropsPhaseType::correctiveMolarVolumeT, "Return the corrective temperature derivative of the molar volume of the phase (in m³/(mol·K)).")
        .def("correctiveMolarVolumeP", &ChemicalPropsPhaseType::correctiveMolarVolumeP, "Return the corrective pressure derivative of the molar volume of the phase (in m³/(mol·Pa)).")
        .def("molarGibbsEnergy", &ChemicalPropsPhaseType::molarGibbsEnergy, "Return the molar Gibbs energy of formation of the phase (in J/mol).")
        .def("molarEnthalpy", &ChemicalPropsPhaseType::molarEnthalpy, "Return the molar enthalpy of formation of the phase (in J/mol).")
        .def("correctiveMolarGibbsEnergy", &ChemicalPropsPhaseType::correctiveMolarGibbsEnergy, "Return the corrective molar Gibbs energy of formation of the phase (in J/mol).")
        .def("correctiveMolarEnthalpy", &ChemicalPropsPhaseType::correctiveMolarEnthalpy, "Return the corrective molar enthalpy of formation of the phase (in J/mol).")
        .def("molarEntropy", &ChemicalPropsPhaseType::molarEntropy, "Return the molar entropy of formation of the phase (in J/(mol·K)).")
        .def("molarInternalEnergy", &ChemicalPropsPhaseType::molarInternalEnergy, "Return the molar internal energy of formation of the phase (in J/mol).")
        .def("molarHelmholtzEnergy", &ChemicalPropsPhaseType::molarHelmholtzEnergy, "Return the molar Helmholtz energy of formation of the phase (in J/mol).")
        .def("molarHeatCapacityConstP", &ChemicalPropsPhaseType::molarHeatCapacityConstP, "Return the molar isobaric heat capacity of the phase (in J/(mol·K)).")
        .def("molarHeatCapacityConstV", &ChemicalPropsPhaseType::molarHeatCapacityConstV, "Return the molar isochoric heat capacity of the phase (in J/(mol·K)).")
        .def("specificVolume", &ChemicalPropsPhaseType::specificVolume, "Return the specific volume of the phase (in m³/kg).")
        .def("specificVolumeT", &ChemicalPropsPhaseType::specificVolumeT, "Return the temperature derivative of the specific volume of the phase (in m³/(kg·K)).")
        .def("specificVolumeP", &ChemicalPropsPhaseType::specificVolumeP, "Return the pressure derivative of the specific volume of the phase (in m³/(kg·Pa)).")
        .def("specificGibbsEnergy", &ChemicalPropsPhaseType::specificGibbsEnergy, "Return the specific Gibbs energy of formation of the phase (in J/kg).")
        .def("specificEnthalpy", &ChemicalPropsPhaseType::specificEnthalpy, "Return the specific enthalpy of formation of the phase (in J/kg).")
        .def("specificEntropy", &ChemicalPropsPhaseType::specificEntropy, "Return the specific entropy of formation of the phase (in J/(kg·K)).")
        .def("specificInternalEnergy", &ChemicalPropsPhaseType::specificInternalEnergy, "Return the specific internal energy of formation of the phase (in J/kg).")
        .def("specificHelmholtzEnergy", &ChemicalPropsPhaseType::specificHelmholtzEnergy, "Return the specific Helmholtz energy of formation of the phase (in J/kg).")
        .def("specificHeatCapacityConstP", &ChemicalPropsPhaseType::specificHeatCapacityConstP, "Return the specific isobaric heat capacity of the phase (in J/(kg·K)).")
        .def("correctiveMolarHeatCapacityConstP", &ChemicalPropsPhaseType::correctiveMolarHeatCapacityConstP, "Return the corrective specific isobaric heat capacity of the phase (in J/(kg·K)).")
        .def("specificHeatCapacityConstV", &ChemicalPropsPhaseType::specificHeatCapacityConstV, "Return the specific isochoric heat capacity of the phase (in J/(kg·K)).")
        .def("density", &ChemicalPropsPhaseType::density, "Return the density of the phase (in kg/m³).")
        .def("amount", &ChemicalPropsPhaseType::amount, "Return the sum of species amounts in the phase (in mol).")
        .def("mass", &ChemicalPropsPhaseType::mass, "Return the sum of species masses in the phase (in kg).")
        .def("volume", &ChemicalPropsPhaseType::volume, "Return the volume of the phase (in m³).")
        .def("volumeT", &ChemicalPropsPhaseType::volumeT, "Return the temperature derivative of the volume of the phase (in m³/K).")
        .def("volumeP", &ChemicalPropsPhaseType::volumeP, "Return the pressure derivative of the volume of the phase (in m³/Pa).")
        .def("gibbsEnergy", &ChemicalPropsPhaseType::gibbsEnergy, "Return the Gibbs energy of the phase (in J).")
        .def("enthalpy", &ChemicalPropsPhaseType::enthalpy, "Return the enthalpy of the phase (in J).")
        .def("entropy", &ChemicalPropsPhaseType::entropy, "Return the entropy of the phase (in J/K).")
        .def("internalEnergy", &ChemicalPropsPhaseType::internalEnergy, "Return the internal energy of the phase (in J).")
        .def("helmholtzEnergy", &ChemicalPropsPhaseType::helmholtzEnergy, "Return the Helmholtz energy of the phase (in J).")
        .def("heatCapacityConstP", &ChemicalPropsPhaseType::heatCapacityConstP, "Return the isobaric heat capacity of the phase (in J/K).")
        .def("heatCapacityConstV", &ChemicalPropsPhaseType::heatCapacityConstV, "Return the isochoric heat capacity of the phase (in J/K).")
        .def("soundSpeed", &ChemicalPropsPhaseType::soundSpeed, "Return the speed of sound in the phase (in m/s).")
        ;
}

void exportChemicalPropsPhase(py::module& m)
{
    createTemplateClassForChemicalPropsPhaseType<ChemicalPropsPhase>(m, "ChemicalPropsPhase")
        .def(py::init<const Phase&>())
        .def("update", &ChemicalPropsPhase::update, "Update the chemical properties of the phase.")
        .def("updateIdeal", &ChemicalPropsPhase::updateIdeal, "Update the chemical properties of the phase using ideal activity models.")
        ;

    createTemplateClassForChemicalPropsPhaseType<ChemicalPropsPhaseRef>(m, "ChemicalPropsPhaseRef")
        .def("update", &ChemicalPropsPhaseRef::update, "Update the chemical properties of the phase.")
        .def("updateIdeal", &ChemicalPropsPhaseRef::updateIdeal, "Update the chemical properties of the phase using ideal activity models.")
        ;

    createTemplateClassForChemicalPropsPhaseType<ChemicalPropsPhaseConstRef>(m, "ChemicalPropsPhaseConstRef")
        ;
}
