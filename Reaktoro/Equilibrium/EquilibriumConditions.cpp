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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace {

/// Throw an error if an input variable has not been registered in the equilibrium specifications.
/// @param inputs The list of registered *w* input variables in the equilibrium specifications (e.g., {"T", "P", "pH"}).
/// @param wid The id of the input variable that needs to be checked in the equilibrium specifications (e.g., "T").
/// @param propertymsg The message about the property being constrained (e.g., "temperature").
auto throwErrorIfNotRegisteredInput(Strings const& inputs, String const& wid, String const& propertymsg) -> void
{
    const auto registered = contains(inputs, wid);
    errorif(!registered, "Cannot set ", propertymsg, " for the equilibrium calculation "
        "because it is not a registered input variable in the equilibrium specifications.");
}

} // namespace

EquilibriumConditions::EquilibriumConditions(ChemicalSystem const& system)
: EquilibriumConditions(EquilibriumSpecs::TP(system))
{}

EquilibriumConditions::EquilibriumConditions(EquilibriumSpecs const& specs)
: msystem(specs.system()),
  C(specs.assembleConservationMatrix()),
  wvars(specs.inputs()),
  pvars(specs.namesControlVariablesP()),
  itemperature_w(specs.indexTemperatureAmongInputVariables()),
  itemperature_p(specs.indexTemperatureAmongControlVariablesP()),
  ipressure_w(specs.indexPressureAmongInputVariables()),
  ipressure_p(specs.indexPressureAmongControlVariablesP())
{
    // Initialize the values of the *w* input variables to NaN to check later if initialized.
    w = constants(specs.numInputs(), NaN);

    // Initialize the default lower and upper bounds for the *p* control variables (-inf and inf respectively).
    plower.setConstant(pvars.size(), -inf);
    pupper.setConstant(pvars.size(),  inf);
}

//=================================================================================================
//
// METHODS TO SPECIFY THERMODYNAMIC CONDITIONS
//
//=================================================================================================

auto EquilibriumConditions::temperature(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "T", "temperature");
    const auto idx = index(wvars, "T");
    w[idx] = units::convert(value, unit, "K");
}

auto EquilibriumConditions::pressure(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "P", "pressure");
    const auto idx = index(wvars, "P");
    w[idx] = units::convert(value, unit, "Pa");
}

auto EquilibriumConditions::volume(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "V", "volume");
    const auto idx = index(wvars, "V");
    w[idx] = units::convert(value, unit, "m3");
}

auto EquilibriumConditions::internalEnergy(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "U", "internal energy");
    const auto idx = index(wvars, "U");
    w[idx] = units::convert(value, unit, "J");
}

auto EquilibriumConditions::enthalpy(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "H", "enthalpy");
    const auto idx = index(wvars, "H");
    w[idx] = units::convert(value, unit, "J");
}

auto EquilibriumConditions::gibbsEnergy(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "G", "Gibbs energy");
    const auto idx = index(wvars, "G");
    w[idx] = units::convert(value, unit, "J");
}

auto EquilibriumConditions::helmholtzEnergy(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "A", "Helmholtz energy");
    const auto idx = index(wvars, "A");
    w[idx] = units::convert(value, unit, "J");
}

auto EquilibriumConditions::entropy(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "S", "entropy");
    const auto idx = index(wvars, "S");
    w[idx] = units::convert(value, unit, "J/K");
}

auto EquilibriumConditions::charge(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "charge", "charge");
    const auto idx = index(wvars, "charge");
    w[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::elementAmount(StringOrIndex const& element, real const& value, String const& unit) -> void
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto elementsymbol = msystem.element(ielement).symbol();
    const auto id = "elementAmount[" + elementsymbol + "]";
    const auto errormsg = "element amount of " + elementsymbol;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::elementAmountInPhase(StringOrIndex const& element, StringOrIndex const& phase, real const& value, String const& unit) -> void
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto elementsymbol = msystem.element(ielement).symbol();
    const auto phasename = msystem.phase(iphase).name();
    const auto id = "elementAmountInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto errormsg = "element amount of " + elementsymbol + " in phase " + phasename;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::elementMass(StringOrIndex const& element, real const& value, String const& unit) -> void
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto elementsymbol = msystem.element(ielement).symbol();
    const auto id = "elementMass[" + elementsymbol + "]";
    const auto errormsg = "element mass of " + elementsymbol;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "kg");
}

auto EquilibriumConditions::elementMassInPhase(StringOrIndex const& element, StringOrIndex const& phase, real const& value, String const& unit) -> void
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto elementsymbol = msystem.element(ielement).symbol();
    const auto phasename = msystem.phase(iphase).name();
    const auto id = "elementMassInPhase[" + elementsymbol + "][" + phasename + "]";
    const auto errormsg = "element mass of " + elementsymbol + " in phase " + phasename;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "kg");
}

auto EquilibriumConditions::phaseAmount(StringOrIndex const& phase, real const& value, String const& unit) -> void
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto phasename = msystem.phase(iphase).name();
    const auto id = "phaseAmount[" + phasename + "]";
    const auto errormsg = "phase amount of " + phasename;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::phaseMass(StringOrIndex const& phase, real const& value, String const& unit) -> void
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto phasename = msystem.phase(iphase).name();
    const auto id = "phaseMass[" + phasename + "]";
    const auto errormsg = "phase mass of " + phasename;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "kg");
}

auto EquilibriumConditions::phaseVolume(StringOrIndex const& phase, real const& value, String const& unit) -> void
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto phasename = msystem.phase(iphase).name();
    const auto id = "phaseVolume[" + phasename + "]";
    const auto errormsg = "phase volume of " + phasename;
    throwErrorIfNotRegisteredInput(wvars, id, errormsg);
    const auto idx = index(wvars, id);
    w[idx] = units::convert(value, unit, "m3");
}

//=================================================================================================
//
// METHODS TO SPECIFY CHEMICAL POTENTIAL CONDITIONS
//
//=================================================================================================

auto EquilibriumConditions::chemicalPotential(String const& substance, real const& value, String const& unit) -> void
{
    const auto pid = "u[" + substance + "]";
    throwErrorIfNotRegisteredInput(wvars, pid, "the chemical potential of " + substance);
    const auto idx = index(wvars, pid);
    w[idx] = units::convert(value, unit, "J/mol");
}

auto EquilibriumConditions::lnActivity(String const& species, real const& value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(wvars, pid, "the activity of " + species);
    const auto idx = index(wvars, pid);
    w[idx] = value;
}

auto EquilibriumConditions::lgActivity(String const& species, real const& value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(wvars, pid, "the activity of " + species);
    const auto idx = index(wvars, pid);
    w[idx] = value * ln10;
}

auto EquilibriumConditions::activity(String const& species, real const& value) -> void
{
    const auto pid = "ln(a[" + species + "])";
    throwErrorIfNotRegisteredInput(wvars, pid, "the activity of " + species);
    const auto idx = index(wvars, pid);
    w[idx] = log(value);
}

auto EquilibriumConditions::fugacity(String const& gas, real const& value, String const& unit) -> void
{
    const auto pid = "f[" + gas + "]";
    throwErrorIfNotRegisteredInput(wvars, pid, "the fugacity of " + gas);
    const auto idx = index(wvars, pid);
    w[idx] = units::convert(value, unit, "bar");
}

auto EquilibriumConditions::pH(real const& value) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "pH", "pH");
    const auto idx = index(wvars, "pH");
    w[idx] = value;
}

auto EquilibriumConditions::pMg(real const& value) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "pMg", "pMg");
    const auto idx = index(wvars, "pMg");
    w[idx] = value;
}

auto EquilibriumConditions::pE(real const& value) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "pE", "pE");
    const auto idx = index(wvars, "pE");
    w[idx] = value;
}

auto EquilibriumConditions::Eh(real const& value, String const& unit) -> void
{
    throwErrorIfNotRegisteredInput(wvars, "Eh", "Eh");
    const auto idx = index(wvars, "Eh");
    w[idx] = value;
}

//=================================================================================================
//
// METHODS FOR SETTING AND GETTING LOWER AND UPPER BOUNDS FOR UNKNOWN VARIABLES
//
//=================================================================================================

auto EquilibriumConditions::setLowerBoundTemperature(double value, String const& unit) -> void
{
    if(itemperature_p < plower.size())
        plower[itemperature_p] = units::convert(value, unit, "K");
}

auto EquilibriumConditions::setUpperBoundTemperature(double value, String const& unit) -> void
{
    if(itemperature_p < pupper.size())
        pupper[itemperature_p] = units::convert(value, unit, "K");
}

auto EquilibriumConditions::setLowerBoundPressure(double value, String const& unit) -> void
{
    if(ipressure_p < plower.size())
        plower[ipressure_p] = units::convert(value, unit, "Pa");
}

auto EquilibriumConditions::setUpperBoundPressure(double value, String const& unit) -> void
{
    if(ipressure_p < pupper.size())
        pupper[ipressure_p] = units::convert(value, unit, "Pa");
}

auto EquilibriumConditions::setLowerBoundTitrant(String const& substance, double value, String const& unit) -> void
{
    const auto idx = index(pvars, substance);
    const auto size = pvars.size();
    errorif(idx >= size, "EquilibriumConditions::setLowerBoundTitrant requires a substance name that was specified in a call to EquilibriumSpecs::openTo.");
    plower[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::setUpperBoundTitrant(String const& substance, double value, String const& unit) -> void
{
    const auto idx = index(pvars, substance);
    const auto size = pvars.size();
    errorif(idx >= size, "EquilibriumConditions::setUpperBoundTitrant requires a substance name that was specified in a call to EquilibriumSpecs::openTo.");
    pupper[idx] = units::convert(value, unit, "mol");
}

auto EquilibriumConditions::setLowerBoundsControlVariablesP(ArrayXdConstRef const& values) -> void
{
    errorif(values.size() != plower.size(), "Expecting in EquilibriumConditions::setLowerBoundsControlVariablesP a vector with same size as that of number of p control variables, ", plower.size(), ", but got instead a vector with size ", values.size(), ".");
    plower = values;
}

auto EquilibriumConditions::setUpperBoundsControlVariablesP(ArrayXdConstRef const& values) -> void
{
    errorif(values.size() != pupper.size(), "Expecting in EquilibriumConditions::setUpperBoundsControlVariablesP a vector with same size as that of number of p control variables, ", pupper.size(), ", but got instead a vector with size ", values.size(), ".");
    pupper = values;
}

auto EquilibriumConditions::lowerBoundsControlVariablesP() const -> ArrayXdConstRef
{
    return plower;
}

auto EquilibriumConditions::upperBoundsControlVariablesP() const -> ArrayXdConstRef
{
    return pupper;
}

//=================================================================================================
//
// METHODS FOR SETTING AND GETTING INPUT VARIABLES
//
//=================================================================================================

auto EquilibriumConditions::set(String const& name, real const& value) -> void
{
    setInputVariable(name, value);
}

auto EquilibriumConditions::setInputVariable(String const& name, real const& value) -> void
{
    const auto idx = index(wvars, name);
    const auto size = wvars.size();
    errorif(idx >= size, "There is no input variable with name `", name, "` in this EquilibriumConditions object.");
    w[idx] = value;
}

auto EquilibriumConditions::setInputVariable(Index index, real const& value) -> void
{
    const auto size = wvars.size();
    errorif(index >= size, "There is no input variable with index ", index, " in this EquilibriumConditions object.");
    w[index] = value;
}

auto EquilibriumConditions::setInputVariables(ArrayXrConstRef const& values) -> void
{
    errorif(values.size() != w.size(), "Expecting in EquilibriumConditions::setInputVariables a vector with same size as that of number of *w* input variables, ", w.size(), ", but got instead a vector with size ", values.size(), ".");
    w = values;
}

auto EquilibriumConditions::inputNames() const -> Strings const&
{
    return wvars;
}

auto EquilibriumConditions::inputValues() const -> ArrayXrConstRef
{
    return w;
}

auto EquilibriumConditions::inputValuesGetOrCompute(ChemicalState const& state0) const -> ArrayXr
{
    // The input values with nan replaced by appropriate values whenever possible
    ArrayXr wvals(w);

    // If temperature is input, but current value is nan, fetch it from state0
    if(itemperature_w < w.size() && std::isnan(w[itemperature_w].val()))
        wvals[itemperature_w] = state0.temperature();

    // If pressure is input, but current value is nan, fetch it from state0
    if(ipressure_w < w.size() && std::isnan(w[ipressure_w].val()))
        wvals[ipressure_w] = state0.pressure();

    // Ensure no other input values are left unspecified! Only temperature and pressure can be inferred at the moment.
    for(auto const& [i, wval] : enumerate(wvals))
        errorif(std::isnan(wval.val()), "You have not specified a value for input `", wvars[i], "` in the EquilibriumConditions object.");

    return wvals;
}

auto EquilibriumConditions::inputValue(String const& name) const -> real const&
{
    auto k = index(wvars, name);
    errorifnot(k < wvars.size(), "Your equilibrium problem specifications do not include an input variable named `", name, "`.");
    return w[k];
}

//=================================================================================================
//
// METHODS TO SPECIFY THE INITIAL COMPOSITIONAL STATE OF THE CHEMICAL SYSTEM BEFORE IT REACTS
//
//=================================================================================================

auto EquilibriumConditions::setInitialComponentAmounts(VectorXdConstRef const& c0) -> void
{
    errorif(c0.rows() != C.rows(), "Expecting a vector of initial amounts of conservative components with size ", C.rows(), " but given one has size ", c0.rows(), " instead.");
    this->c0 = c0;
}

auto EquilibriumConditions::setInitialComponentAmountsFromSpeciesAmounts(VectorXdConstRef const& n0) -> void
{
    errorif(n0.rows() != C.cols(), "Expecting a vector of initial amounts of species with size ", C.cols(), " but given one has size ", n0.rows(), " instead.");
    c0 = C * n0;
}

auto EquilibriumConditions::setInitialComponentAmountsFromState(ChemicalState const& state0) -> void
{
    const VectorXdConstRef n0 = state0.speciesAmounts();
    c0 = C * n0;
}

auto EquilibriumConditions::initialComponentAmounts() const -> ArrayXdConstRef
{
    return c0;
}

auto EquilibriumConditions::initialComponentAmountsGetOrCompute(VectorXdConstRef const& n0) const -> ArrayXd
{
    errorif(n0.rows() != C.cols(), "Expecting a vector of initial amounts of species with size ", C.cols(), " but given one has size ", n0.rows(), " instead.");
    return c0.rows() != 0 ? c0 : ArrayXd(C * n0);
}

auto EquilibriumConditions::initialComponentAmountsGetOrCompute(ChemicalState const& state0) const -> ArrayXd
{
    const VectorXdConstRef n0 = state0.speciesAmounts();
    return c0.rows() != 0 ? c0 : ArrayXd(C * n0);
}

//=================================================================================================
//
// MISCELLANEOUS METHODS
//
//=================================================================================================

auto EquilibriumConditions::system() const -> ChemicalSystem const&
{
    return msystem;
}

} // namespace Reaktoro
