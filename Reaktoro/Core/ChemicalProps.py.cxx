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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
using namespace Reaktoro;

void exportChemicalProps(py::module& m)
{
    py::class_<ChemicalProps>(m, "ChemicalProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def("update", py::overload_cast<const ChemicalState&>(&ChemicalProps::update), "Update the chemical properties of the system.")
        .def("update", py::overload_cast<const real&, const real&, ArrayXrConstRef>(&ChemicalProps::update), "Update the chemical properties of the system.")
        .def("update", py::overload_cast<ArrayXrConstRef>(&ChemicalProps::update), "Update the chemical properties of the system with serialized data.")
        .def("update", py::overload_cast<ArrayXdConstRef>(&ChemicalProps::update), "Update the chemical properties of the system with serialized data.")
        .def("updateIdeal", py::overload_cast<const ChemicalState&>(&ChemicalProps::updateIdeal), "Update the chemical properties of the system using ideal activity models.")
        .def("updateIdeal", py::overload_cast<const real&, const real&, ArrayXrConstRef>(&ChemicalProps::updateIdeal), "Update the chemical properties of the system using ideal activity models.")
        .def("system", &ChemicalProps::system, return_internal_ref, "Return the chemical system associated with these chemical properties.")
        .def("phaseProps", py::overload_cast<StringOrIndex>(&ChemicalProps::phaseProps), return_internal_ref, "Return the chemical properties of a phase with given index.")
        .def("phaseProps", py::overload_cast<StringOrIndex>(&ChemicalProps::phaseProps, py::const_), return_internal_ref, "Return the chemical properties of a phase with given index.")
        .def("temperature", &ChemicalProps::temperature, "Return the temperature of the system (in K).")
        .def("pressure", &ChemicalProps::pressure, "Return the pressure of the system (in Pa).")
        .def("charge", &ChemicalProps::charge, "Return the amount of electric charge in the system (in mol).")
        .def("chargeInPhase", &ChemicalProps::chargeInPhase, "Return the amount of electric charge in a phase of the system (in mol).")
        .def("chargeAmongSpecies", &ChemicalProps::chargeAmongSpecies, "Return the amount of electric charge among a group of species in the system (in mol).")
        .def("elementAmount", &ChemicalProps::elementAmount, "Return the amount of an element in the system (in mol).")
        .def("elementAmountInPhase", &ChemicalProps::elementAmountInPhase, "Return the amount of an element in the system (in mol).")
        .def("elementAmountAmongSpecies", &ChemicalProps::elementAmountAmongSpecies, "Return the amount of an element among a group of species in the system (in mol).")
        .def("elementMass", &ChemicalProps::elementMass, "Return the mass of an element in the system (in kg).")
        .def("elementMassInPhase", &ChemicalProps::elementMassInPhase, "Return the mass of an element in the system (in kg).")
        .def("elementMassAmongSpecies", &ChemicalProps::elementMassAmongSpecies, "Return the mass of an element among a group of species in the system (in kg).")
        .def("elementAmounts", &ChemicalProps::elementAmounts, "Return the amounts of the elements in the system (in mol).")
        .def("elementAmountsInPhase", &ChemicalProps::elementAmountsInPhase, "Return the amounts of the elements in a phase of the system (in mol).")
        .def("elementAmountsAmongSpecies", &ChemicalProps::elementAmountsAmongSpecies, "Return the amounts of the elements among a group of species in the system (in mol).")
        .def("speciesAmount", &ChemicalProps::speciesAmount, "Return the amount of a species in the system (in mol).")
        .def("speciesMass", &ChemicalProps::speciesMass, "Return the mass of a species in the system.")
        .def("speciesMoleFraction", &ChemicalProps::speciesMoleFraction, "Return the mole fraction of a species in the system.")
        .def("speciesConcentration", &ChemicalProps::speciesConcentration, "Return the concentration (activity divided by activity coefficient) of a species in the system.")
        .def("speciesConcentrationLg", &ChemicalProps::speciesConcentrationLg, "Return the lg concentration (activity divided by activity coefficient) of a species in the system.")
        .def("speciesConcentrationLn", &ChemicalProps::speciesConcentrationLn, "Return the ln concentration (activity divided by activity coefficient) of a species in the system.")
        .def("speciesActivityCoefficient", &ChemicalProps::speciesActivityCoefficient, "Return the activity coefficient of a species in the system.")
        .def("speciesActivityCoefficientLg", &ChemicalProps::speciesActivityCoefficientLg, "Return the lg activity coefficient of a species in the system.")
        .def("speciesActivityCoefficientLn", &ChemicalProps::speciesActivityCoefficientLn, "Return the ln activity coefficient of a species in the system.")
        .def("speciesActivity", &ChemicalProps::speciesActivity, "Return the activity of a species in the system.")
        .def("speciesActivityLg", &ChemicalProps::speciesActivityLg, "Return the lg activity of a species in the system.")
        .def("speciesActivityLn", &ChemicalProps::speciesActivityLn, "Return the ln activity of a species in the system.")
        .def("speciesChemicalPotential", &ChemicalProps::speciesChemicalPotential, "Return the chemical potential of a species in the system.")
        .def("speciesStandardVolume", &ChemicalProps::speciesStandardVolume, "Return the standard partial molar volume of a species in the system (in m³/mol).")
        .def("speciesStandardVolumeT", &ChemicalProps::speciesStandardVolumeT, "Return the temperature derivative of the standard partial molar volume of a species in the system (in m³/(mol·K)).")
        .def("speciesStandardVolumeP", &ChemicalProps::speciesStandardVolumeP, "Return the pressure derivative of the standard partial molar volume of a species in the system (in m³/(mol·Pa)).")
        .def("speciesStandardGibbsEnergy", &ChemicalProps::speciesStandardGibbsEnergy, "Return the standard partial molar Gibbs energy of formation of a species in the system (in J/mol).")
        .def("speciesStandardEnthalpy", &ChemicalProps::speciesStandardEnthalpy, "Return the standard partial molar enthalpy of formation of a species in the system (in J/mol).")
        .def("speciesStandardEntropy", &ChemicalProps::speciesStandardEntropy, "Return the standard partial molar entropy of formation of the species a the system (in J/(mol·K)).")
        .def("speciesStandardInternalEnergy", &ChemicalProps::speciesStandardInternalEnergy, "Return the standard partial molar internal energy of formation of a species in the system (in J/mol).")
        .def("speciesStandardHelmholtzEnergy", &ChemicalProps::speciesStandardHelmholtzEnergy, "Return the standard partial molar Helmholtz energy of formation of a species in the system (in J/mol).")
        .def("speciesStandardHeatCapacityConstP", &ChemicalProps::speciesStandardHeatCapacityConstP, "Return the standard partial molar isobaric heat capacity of the species a the system (in J/(mol·K)).")
        .def("speciesStandardHeatCapacityConstV", &ChemicalProps::speciesStandardHeatCapacityConstV, "Return the standard partial molar isochoric heat capacity of the species a the system (in J/(mol·K)).")
        .def("speciesAmounts", &ChemicalProps::speciesAmounts, return_internal_ref, "Return the amounts of the species in the system (in mol).")
        .def("speciesMasses", &ChemicalProps::speciesMasses, "Return the masses of the species in the system (in kg).")
        .def("speciesMoleFractions", &ChemicalProps::speciesMoleFractions, return_internal_ref, "Return the mole fractions of the species in the system.")
        .def("speciesConcentrationsLn", &ChemicalProps::speciesConcentrationsLn, "Return the ln concentrations (activity divided by activity coefficient) of the species in the system.")
        .def("speciesActivityCoefficientsLn", &ChemicalProps::speciesActivityCoefficientsLn, return_internal_ref, "Return the ln activity coefficients of the species in the system.")
        .def("speciesActivitiesLn", &ChemicalProps::speciesActivitiesLn, return_internal_ref, "Return the ln activities of the species in the system.")
        .def("speciesChemicalPotentials", &ChemicalProps::speciesChemicalPotentials, return_internal_ref, "Return the chemical potentials of the species in the system (in J/mol).")
        .def("speciesStandardVolumes", &ChemicalProps::speciesStandardVolumes, return_internal_ref, "Return the standard partial molar volumes of the species in the system (in m³/mol).")
        .def("speciesStandardVolumesT", &ChemicalProps::speciesStandardVolumesT, return_internal_ref, "Return the temperature derivative of the standard molar volumes of the species in the system (in m³/(mol·K)).")
        .def("speciesStandardVolumesP", &ChemicalProps::speciesStandardVolumesP, return_internal_ref, "Return the pressure derivative of the standard molar volumes of the species in the system (in m³/(mol·Pa)).")
        .def("speciesStandardGibbsEnergies", &ChemicalProps::speciesStandardGibbsEnergies, return_internal_ref, "Return the standard partial molar Gibbs energies of formation of the species in the system (in J/mol).")
        .def("speciesStandardEnthalpies", &ChemicalProps::speciesStandardEnthalpies, return_internal_ref, "Return the standard partial molar enthalpies of formation of the species in the system (in J/mol).")
        .def("speciesStandardEntropies", &ChemicalProps::speciesStandardEntropies, "Return the standard partial molar entropies of formation of the species in the system (in J/(mol·K)).")
        .def("speciesStandardInternalEnergies", &ChemicalProps::speciesStandardInternalEnergies, "Return the standard partial molar internal energies of formation of the species in the system (in J/mol).")
        .def("speciesStandardHelmholtzEnergies", &ChemicalProps::speciesStandardHelmholtzEnergies, "Return the standard partial molar Helmholtz energies of formation of the species in the system (in J/mol).")
        .def("speciesStandardHeatCapacitiesConstP", &ChemicalProps::speciesStandardHeatCapacitiesConstP, return_internal_ref, "Return the standard partial molar isobaric heat capacities of the species in the system (in J/(mol·K)).")
        .def("speciesStandardHeatCapacitiesConstV", &ChemicalProps::speciesStandardHeatCapacitiesConstV, return_internal_ref, "Return the standard partial molar isochoric heat capacities of the species in the system (in J/(mol·K)).")
        .def("molarVolume", &ChemicalProps::molarVolume, "Return the molar volume of the system (in m³/mol).")
        .def("molarVolumeT", &ChemicalProps::molarVolumeT, "Return the temperature derivative of the molar volume of the system (in m³/(mol·K)).")
        .def("molarVolumeP", &ChemicalProps::molarVolumeP, "Return the pressure derivative of the molar volume of the system (in m³/(mol·Pa)).")
        .def("molarGibbsEnergy", &ChemicalProps::molarGibbsEnergy, "Return the molar Gibbs energy of formation of the system (in J/mol).")
        .def("molarEnthalpy", &ChemicalProps::molarEnthalpy, "Return the molar enthalpy of formation of the system (in J/mol).")
        .def("molarEntropy", &ChemicalProps::molarEntropy, "Return the molar entropy of formation of the system (in J/(mol·K)).")
        .def("molarInternalEnergy", &ChemicalProps::molarInternalEnergy, "Return the molar internal energy of formation of the system (in J/mol).")
        .def("molarHelmholtzEnergy", &ChemicalProps::molarHelmholtzEnergy, "Return the molar Helmholtz energy of formation of the system (in J/mol).")
        .def("molarHeatCapacityConstP", &ChemicalProps::molarHeatCapacityConstP, "Return the molar isobaric heat capacity of the system (in J/(mol·K)).")
        .def("molarHeatCapacityConstV", &ChemicalProps::molarHeatCapacityConstV, "Return the molar isochoric heat capacity of the system (in J/(mol·K)).")
        .def("specificVolume", &ChemicalProps::specificVolume, "Return the specific volume of the system (in m³/kg).")
        .def("specificVolumeT", &ChemicalProps::specificVolumeT, "Return the temperature derivative of the specific volume of the system (in m³/(kg·K)).")
        .def("specificVolumeP", &ChemicalProps::specificVolumeP, "Return the pressure derivative of the specific volume of the system (in m³/(kg·Pa)).")
        .def("specificGibbsEnergy", &ChemicalProps::specificGibbsEnergy, "Return the specific Gibbs energy of formation of the system (in J/kg).")
        .def("specificEnthalpy", &ChemicalProps::specificEnthalpy, "Return the specific enthalpy of formation of the system (in J/kg).")
        .def("specificEntropy", &ChemicalProps::specificEntropy, "Return the specific entropy of formation of the system (in J/(kg·K)).")
        .def("specificInternalEnergy", &ChemicalProps::specificInternalEnergy, "Return the specific internal energy of formation of the system (in J/kg).")
        .def("specificHelmholtzEnergy", &ChemicalProps::specificHelmholtzEnergy, "Return the specific Helmholtz energy of formation of the system (in J/kg).")
        .def("specificHeatCapacityConstP", &ChemicalProps::specificHeatCapacityConstP, "Return the specific isobaric heat capacity of the system (in J/(kg·K)).")
        .def("specificHeatCapacityConstV", &ChemicalProps::specificHeatCapacityConstV, "Return the specific isochoric heat capacity of the system (in J/(kg·K)).")
        .def("density", &ChemicalProps::density, "Return the density of the system (in kg/m³).")
        .def("amount", &ChemicalProps::amount, "Return the sum of species amounts in the system (in mol).")
        .def("mass", &ChemicalProps::mass, "Return the sum of species masses in the system (in kg).")
        .def("volume", &ChemicalProps::volume, "Return the volume of the system (in m³).")
        .def("volumeT", &ChemicalProps::volumeT, "Return the temperature derivative of the volume of the system (in m³/K).")
        .def("volumeP", &ChemicalProps::volumeP, "Return the pressure derivative of the volume of the system (in m³/Pa).")
        .def("gibbsEnergy", &ChemicalProps::gibbsEnergy, "Return the Gibbs energy of formation of the system (in J).")
        .def("enthalpy", &ChemicalProps::enthalpy, "Return the enthalpy of formation of the system (in J).")
        .def("entropy", &ChemicalProps::entropy, "Return the entropy of formation of the system (in J/K).")
        .def("internalEnergy", &ChemicalProps::internalEnergy, "Return the internal energy of formation of the system (in J).")
        .def("helmholtzEnergy", &ChemicalProps::helmholtzEnergy, "Return the Helmholtz energy of formation of the system (in J).")
        .def("heatCapacityConstP", &ChemicalProps::heatCapacityConstP, "Return the isobaric heat capacity of the system (in J/K).")
        .def("heatCapacityConstV", &ChemicalProps::heatCapacityConstV, "Return the isochoric heat capacity of the system (in J/K).")
        .def("indicesPhasesWithState", &ChemicalProps::indicesPhasesWithState, "Return the indices of the phases in a given state of matter.")
        .def("indicesPhasesWithStates", &ChemicalProps::indicesPhasesWithStates, "Return the indices of the phases with one of the given states of matter.")
        .def("indicesPhasesWithFluidState", &ChemicalProps::indicesPhasesWithFluidState, "Return the indices of the phases in liquid, gaseous, or supercritical states.")
        .def("indicesPhasesWithSolidState", &ChemicalProps::indicesPhasesWithSolidState, "Return the indices of the phases in solid states.")
        .def("indicesSpeciesInPhasesWithState", &ChemicalProps::indicesSpeciesInPhasesWithState, "Return the indices of the phases with a given state of matter.")
        .def("indicesSpeciesInPhasesWithStates", &ChemicalProps::indicesSpeciesInPhasesWithStates, "Return the indices of the phases with one of the given states of matter.")
        .def("indicesSpeciesInPhasesWithFluidState", &ChemicalProps::indicesSpeciesInPhasesWithFluidState, "Return the indices of the phases with liquid, gaseous, or supercritical states.")
        .def("indicesSpeciesInPhasesWithSolidState", &ChemicalProps::indicesSpeciesInPhasesWithSolidState, "Return the indices of the phases with solid states.")
        .def("output", py::overload_cast<std::ostream&>(&ChemicalProps::output, py::const_), "Output the chemical properties of the system to a stream.")
        .def("output", py::overload_cast<const String&>(&ChemicalProps::output, py::const_), "Output the chemical properties of the system to a file.")
        .def("__repr__", [](const ChemicalProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
