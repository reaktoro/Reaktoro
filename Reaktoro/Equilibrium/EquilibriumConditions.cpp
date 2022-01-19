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

#include "EquilibriumConditions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

/// Throw an error if an input variable has not been registered in the equilibrium specifications.
/// @param inputs The list of registered input variables in the equilibrium specifications (e.g., {"T", "P", "pH"}).
/// @param wid The id of the input variable that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredInput(const Strings& inputs, const String& wid, const String& propertymsg) -> void
{
    const auto registered = contains(inputs, wid);
    errorif(!registered, "Cannot set ", propertymsg, " for the equilibrium calculation "
        "because it is not a registered input variable in the equilibrium specifications.");
}

} // namespace

EquilibriumConditions::EquilibriumConditions(const EquilibriumSpecs& specs)
: m_system(specs.system()),
  m_inputs(specs.inputs()),
  m_control_variables_p(specs.namesControlVariablesP()),
  m_idxT(specs.indexControlVariableTemperature()),
  m_idxP(specs.indexControlVariablePressure())
{
    // Initialize the values of the input variables to zero
    m_inputs_values = zeros(specs.numInputs());

    // Initialize the values of the input variables that are model parameters to their current values
    for(auto [i, idx] : enumerate(specs.indicesParams()))
        m_inputs_values[idx] = specs.params()[i].value();

    // Initialize the default lower and upper bounds for the *p* control variables (-inf and inf respectively).
    m_plower.setConstant(m_control_variables_p.size(), -inf);
    m_pupper.setConstant(m_control_variables_p.size(),  inf);
}

auto EquilibriumConditions::temperature(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "T", "temperature");
    value = units::convert(value, unit, "K");
    const auto idx = index(m_inputs, "T");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pressure(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "P", "pressure");
    value = units::convert(value, unit, "Pa");
    const auto idx = index(m_inputs, "P");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::volume(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "V", "volume");
    value = units::convert(value, unit, "m3");
    const auto idx = index(m_inputs, "V");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::internalEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "U", "internal energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "U");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::enthalpy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "H", "enthalpy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "H");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::gibbsEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "G", "Gibbs energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "G");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::helmholtzEnergy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "A", "Helmholtz energy");
    value = units::convert(value, unit, "J");
    const auto idx = index(m_inputs, "A");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::entropy(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "S", "entropy");
    value = units::convert(value, unit, "J/K");
    const auto idx = index(m_inputs, "S");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::charge(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "charge", "charge");
    value = units::convert(value, unit, "mol");
    const auto idx = index(m_inputs, "charge");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::elementAmount(const StringOrIndex& element, real value, String unit) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto id = "elementAmount[" + elementsymbol + "]";
    const auto errormsg = "element amount of " + elementsymbol;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "mol");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::elementAmountInPhase(const StringOrIndex& element, const StringOrIndex& phase, real value, String unit) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto phasename = m_system.phase(iphase).name();
    const auto id = "elementAmountInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto errormsg = "element amount of " + elementsymbol + " in phase " + phasename;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "mol");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::elementMass(const StringOrIndex& element, real value, String unit) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto id = "elementMass[" + elementsymbol + "]";
    const auto errormsg = "element mass of " + elementsymbol;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "kg");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::elementMassInPhase(const StringOrIndex& element, const StringOrIndex& phase, real value, String unit) -> void
{
    const auto ielement = detail::resolveElementIndex(m_system, element);
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto elementsymbol = m_system.element(ielement).symbol();
    const auto phasename = m_system.phase(iphase).name();
    const auto id = "elementMassInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto errormsg = "element mass of " + elementsymbol + " in phase " + phasename;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "kg");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::phaseAmount(const StringOrIndex& phase, real value, String unit) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    const auto id = "phaseAmount[" + phasename + "]";
    const auto errormsg = "phase amount of " + phasename;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "mol");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::phaseMass(const StringOrIndex& phase, real value, String unit) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    const auto id = "phaseMass[" + phasename + "]";
    const auto errormsg = "phase mass of " + phasename;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "kg");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::phaseVolume(const StringOrIndex& phase, real value, String unit) -> void
{
    const auto iphase = detail::resolvePhaseIndex(m_system, phase);
    const auto phasename = m_system.phase(iphase).name();
    const auto id = "phaseVolume[" + phasename + "]";
    const auto errormsg = "phase volume of " + phasename;
    throwErrorIfNotRegisteredInput(m_inputs, id, errormsg);
    value = units::convert(value, unit, "m3");
    const auto idx = index(m_inputs, id);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::chemicalPotential(String substance, real value, String unit) -> void
{
    const auto pid = "u[" + substance + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the chemical potential of " + substance);
    value = units::convert(value, unit, "J/mol");
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::lnActivity(String species, real value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::lgActivity(String species, real value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value * ln10;
}

auto EquilibriumConditions::activity(String species, real value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the activity of " + species);
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = log(value);
}

auto EquilibriumConditions::fugacity(String gas, real value, String unit) -> void
{
    const auto pid = "f[" + gas + "]";
    throwErrorIfNotRegisteredInput(m_inputs, pid, "the fugacity of " + gas);
    value = units::convert(value, unit, "bar");
    const auto idx = index(m_inputs, pid);
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pH(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pH", "pH");
    const auto idx = index(m_inputs, "pH");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pMg(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pMg", "pMg");
    const auto idx = index(m_inputs, "pMg");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::pE(real value) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "pE", "pE");
    const auto idx = index(m_inputs, "pE");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::Eh(real value, String unit) -> void
{
    throwErrorIfNotRegisteredInput(m_inputs, "Eh", "Eh");
    const auto idx = index(m_inputs, "Eh");
    m_inputs_values[idx] = value;
}

auto EquilibriumConditions::setLowerBoundTemperature(double value, String unit) -> void
{
    if(m_idxT < m_plower.size())
        m_plower[m_idxT] = units::convert(value, unit, "K");
}

auto EquilibriumConditions::setUpperBoundTemperature(double value, String unit) -> void
{
    if(m_idxT < m_pupper.size())
        m_pupper[m_idxT] = units::convert(value, unit, "K");
}

auto EquilibriumConditions::setLowerBoundPressure(double value, String unit) -> void
{
    if(m_idxP < m_plower.size())
        m_plower[m_idxP] = units::convert(value, unit, "Pa");
}

auto EquilibriumConditions::setUpperBoundPressure(double value, String unit) -> void
{
    if(m_idxP < m_pupper.size())
        m_pupper[m_idxP] = units::convert(value, unit, "Pa");
}

auto EquilibriumConditions::setLowerBoundTitrant(String substance, double value, String unit) -> void
{
    const auto idx = index(m_control_variables_p, substance);
    const auto size = m_control_variables_p.size();
    errorif(idx >= size, "EquilibriumConditions::setLowerBoundTitrant requires a substance name that was specified in a call to EquilibriumSpecs::openTo.");
    m_plower[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::setUpperBoundTitrant(String substance, double value, String unit) -> void
{
    const auto idx = index(m_control_variables_p, substance);
    const auto size = m_control_variables_p.size();
    errorif(idx >= size, "EquilibriumConditions::setUpperBoundTitrant requires a substance name that was specified in a call to EquilibriumSpecs::openTo.");
    m_pupper[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::set(const String& input, const real& val) -> void
{
    const auto idx = index(m_inputs, input);
    const auto size = m_inputs.size();
    errorif(idx >= size, "There is no input variable with name `", input, "` in this EquilibriumConditions object.");
    m_inputs_values[idx] = val;
}

auto EquilibriumConditions::system() const -> const ChemicalSystem&
{
    return m_system;
}

auto EquilibriumConditions::inputNames() const -> const Strings&
{
    return m_inputs;
}

auto EquilibriumConditions::inputValues() const -> VectorXrConstRef
{
    return m_inputs_values;
}

auto EquilibriumConditions::lowerBoundsControlVariablesP() const -> VectorXdConstRef
{
    return m_plower;
}

auto EquilibriumConditions::upperBoundsControlVariablesP() const -> VectorXdConstRef
{
    return m_pupper;
}

} // namespace Reaktoro
