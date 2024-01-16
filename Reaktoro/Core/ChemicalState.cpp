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

#include "ChemicalState.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Optima includes
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

//=================================================================================================
//
// ChemicalState
//
//=================================================================================================

struct ChemicalState::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The properties related to an equilibrium state.
    Equilibrium equilibrium;

    /// The chemical properties of the system associated to this chemical state.
    ChemicalProps props;

    /// The temperature state of the chemical system (in K)
    real T = 298.15;

    /// The pressure state of the chemical system (in Pa)
    real P = 1.0e+05;

    /// The amounts of the chemical species (in mol)
    ArrayXr n;

    /// Construct a ChemicalState::Impl instance with given chemical system.
    Impl(ChemicalSystem const& system)
    : system(system), equilibrium(system), props(system)
    {
        n.setConstant(system.species().size(), 1e-16); // set small positive value for initial species amounts
    }

    auto temperature(real const& val) -> void
    {
        errorif(val <= 0.0, "Expecting a positive temperature value, but got ", val, " K.");
        T = val;
    }

    auto temperature(real val, Chars unit) -> void
    {
        temperature(units::convert(val, unit, "K"));
    }

    auto pressure(real const& val) -> void
    {
        errorif(val <= 0.0, "Expecting a positive pressure value, but got ", val, " Pa.");
        P = val;
    }

    auto pressure(real val, Chars unit) -> void
    {
        errorif(val <= 0.0, "Expecting a positive pressure value, but got ", val, " ", unit, ".");
        pressure(units::convert(val, unit, "Pa"));
    }

    // --------------------------------------------------------------------------------------------
    // METHODS FOR SETTING THE AMOUNT OR MASS OF SPECIES
    // --------------------------------------------------------------------------------------------

    auto setSpeciesAmounts(real const& val) -> void
    {
        errorif(val < 0.0, "It is not possible to set a negative value for the species amounts.");
        n.fill(val);
    }

    auto setSpeciesAmounts(ArrayXrConstRef const& values) -> void
    {
        errorif(n.size() != values.size(), "Expecting given vector of species amounts to be of size ", n.size(), " but its size is ", values.size(), ".");
        n = values;
    }

    auto setSpeciesAmounts(ArrayXdConstRef const& values) -> void
    {
        errorif(n.size() != values.size(), "Expecting given vector of species amounts to be of size ", n.size(), " but its size is ", values.size(), ".");
        n = values;
    }

    auto setSpeciesAmount(Index ispecies, real const& amount) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " mol.");
        errorif(ispecies >= system.species().size(), "Given species index ", ispecies, " is greater than number of species, `", system.species().size(), ".");
        n[ispecies] = amount;
    }

    auto setSpeciesAmount(StringOrIndex const& species, real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", stringfy(species), "`.");
        n[ispecies] = units::convert(amount, unit, "mol");
    }

    auto setSpeciesMass(StringOrIndex const& species, real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", stringfy(species), "`.");
        n[ispecies] = units::convert(mass, unit, "kg") / system.species(ispecies).molarMass();
    }

    auto set(StringOrIndex const& species, real value, Chars unit) -> void
    {
        errorif(value < 0.0, "Expecting a non-negative amount/mass value, but got ", value, " ", unit);
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        const auto numspecies = system.species().size();
        errorif(ispecies >= numspecies, "Could not find a species in the system with index or name `", stringfy(species), "`.");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] = amount;
    }

    auto add(StringOrIndex const& species, real value, Chars unit) -> void
    {
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        const auto numspecies = system.species().size();
        errorif(ispecies >= numspecies, "Could not find a species in the system with index or name `", stringfy(species), "`.");
        const auto amount = detail::computeSpeciesAmount(system, ispecies, value, unit);
        n[ispecies] += amount;
        errorif(n[ispecies] < 0.0, "It is not possible to add a negative species amount (", value, " ", unit, ") that produces a negative amount for the species.");
    }

    // --------------------------------------------------------------------------------------------
    // METHODS FOR GETTING THE AMOUNT OR MASS OF SPECIES, ELEMENTS, AND CHARGE
    // --------------------------------------------------------------------------------------------

    auto speciesAmountsInPhase(StringOrIndex const& phase) const -> ArrayXrConstRef
    {
        const auto iphase = detail::resolvePhaseIndexOrRaiseError(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        return n.segment(start, size);
    }

    auto speciesAmount(StringOrIndex const& species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", stringfy(species), "`.");
        return n[ispecies];
    }

    auto speciesMass(StringOrIndex const& species) const -> real
    {
        const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(system, species);
        errorif(ispecies >= system.species().size(), "Could not find a species in the system with index or name `", stringfy(species), "`.");
        return n[ispecies] * system.species(ispecies).molarMass();
    }

    auto componentAmounts() const -> ArrayXr
    {
        auto const& A = system.formulaMatrix();
        return A * n.matrix();
    }

    auto elementAmounts() const -> ArrayXr
    {
        auto const& Ae = system.formulaMatrixElements();
        return Ae * n.matrix();
    }

    auto charge() const -> real
    {
        auto const& Az = system.formulaMatrixCharge();
        return (Az * n.matrix())[0];
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE AMOUNTS OF SPECIES IN THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleSpeciesAmounts(double scalar) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        n *= scalar;
    }

    auto scaleSpeciesAmounts(double scalar, Indices const& indices) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        n(indices) *= scalar;
    }

    auto scaleSpeciesAmountsInPhase(StringOrIndex const& phase, double scalar) -> void
    {
        errorif(scalar < 0.0, "Expecting a non-negative scaling factor, but got ", scalar);
        const auto iphase = detail::resolvePhaseIndexOrRaiseError(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", stringfy(phase));
        const auto start = system.phases().numSpeciesUntilPhase(iphase);
        const auto size = system.phase(iphase).species().size();
        n.segment(start, size) *= scalar;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE AMOUNT OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleAmount(real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        amount = units::convert(amount, unit, "mol");
        props.update(T, P, n);
        const auto current_amount = props.amount();
        const auto scalar = (current_amount != 0.0) ? amount/current_amount : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseAmount(StringOrIndex const& phase, real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        amount = units::convert(amount, unit, "mol");
        const auto iphase = detail::resolvePhaseIndexOrRaiseError(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", stringfy(phase));
        props.update(T, P, n);
        const auto current_amount = props.phaseProps(iphase).amount();
        const auto scalar = (current_amount != 0.0) ? amount/current_amount : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidAmount(real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        amount = units::convert(amount, unit, "mol");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_amount =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).amount(); });
        auto const& factor = current_fluid_amount > 0.0 ? amount / current_fluid_amount : real(0.0);
        auto const& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        n(ifluidspecies) *= factor;
    }

    auto scaleSolidAmount(real amount, Chars unit) -> void
    {
        errorif(amount < 0.0, "Expecting a non-negative amount value, but got ", amount, " ", unit);
        amount = units::convert(amount, unit, "mol");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_amount =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).amount(); });
        auto const& factor = current_solid_amount > 0.0 ? amount / current_solid_amount : real(0.0);
        auto const& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        n(isolidspecies) *= factor;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE MASS OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto current_mass = props.mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseMass(StringOrIndex const& phase, real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        const auto iphase = detail::resolvePhaseIndexOrRaiseError(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", stringfy(phase));
        props.update(T, P, n);
        const auto current_mass = props.phaseProps(iphase).mass();
        const auto scalar = (current_mass != 0.0) ? mass/current_mass : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_mass =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        auto const& factor = current_fluid_mass > 0.0 ? mass / current_fluid_mass : real(0.0);
        auto const& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        n(ifluidspecies) *= factor;
    }

    auto scaleSolidMass(real mass, Chars unit) -> void
    {
        errorif(mass < 0.0, "Expecting a non-negative mass value, but got ", mass, " ", unit);
        mass = units::convert(mass, unit, "kg");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_mass =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).mass(); });
        auto const& factor = current_solid_mass > 0.0 ? mass / current_solid_mass : real(0.0);
        auto const& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        n(isolidspecies) *= factor;
    }

    // --------------------------------------------------------------------------------------------
    // METHODS TO SCALE THE VOLUME OF THE SYSTEM OR PART OF IT
    // --------------------------------------------------------------------------------------------

    auto scaleVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto current_volume = props.volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmounts(scalar);
    }

    auto scalePhaseVolume(StringOrIndex const& phase, real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        const auto iphase = detail::resolvePhaseIndexOrRaiseError(system, phase);
        errorif(iphase >= system.phases().size(), "Could not find a phase in the system with index or name `", stringfy(phase));
        props.update(T, P, n);
        const auto current_volume = props.phaseProps(iphase).volume();
        const auto scalar = (current_volume != 0.0) ? volume/current_volume : real(0.0);
        scaleSpeciesAmountsInPhase(iphase, scalar);
    }

    auto scaleFluidVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto ifluidphases = props.indicesPhasesWithFluidState();
        const auto current_fluid_volume =
            Reaktoro::sum(ifluidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        auto const& factor = current_fluid_volume > 0.0 ? volume / current_fluid_volume : real(0.0);
        auto const& ifluidspecies = system.phases().indicesSpeciesInPhases(ifluidphases);
        n(ifluidspecies) *= factor;
    }

    auto scaleSolidVolume(real volume, Chars unit) -> void
    {
        errorif(volume < 0.0, "Expecting a non-negative volume value, but got ", volume, " ", unit);
        volume = units::convert(volume, unit, "m3");
        props.update(T, P, n);
        const auto isolidphases = props.indicesPhasesWithSolidState();
        const auto current_solid_volume =
            Reaktoro::sum(isolidphases, [&](auto i) { return props.phaseProps(i).volume(); });
        auto const& factor = current_solid_volume > 0.0 ? volume / current_solid_volume : real(0.0);
        auto const& isolidspecies = system.phases().indicesSpeciesInPhases(isolidphases);
        n(isolidspecies) *= factor;
    }
};

ChemicalState::ChemicalState(ChemicalSystem const& system)
: pimpl(new Impl(system))
{}

ChemicalState::ChemicalState(ChemicalState const& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalState::~ChemicalState()
{}

auto ChemicalState::operator=(ChemicalState other) -> ChemicalState&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING/GETTING TEMPERATURE
// --------------------------------------------------------------------------------------------

auto ChemicalState::setTemperature(real const& value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::setTemperature(real value, Chars unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::temperature(real const& value) -> void
{
    pimpl->temperature(value);
}

auto ChemicalState::temperature(real value, Chars unit) -> void
{
    pimpl->temperature(value, unit);
}

auto ChemicalState::temperature() const -> real
{
    return pimpl->T;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING/GETTING PRESSURE
// --------------------------------------------------------------------------------------------

auto ChemicalState::setPressure(real const& value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::setPressure(real value, Chars unit) -> void
{
    pimpl->pressure(value, unit);
}

auto ChemicalState::pressure(real const& value) -> void
{
    pimpl->pressure(value);
}

auto ChemicalState::pressure(real value, Chars unit) -> void
{
    pimpl->pressure(value, unit);
}

auto ChemicalState::pressure() const -> real
{
    return pimpl->P;
}

// --------------------------------------------------------------------------------------------
// METHODS FOR SETTING THE AMOUNT OR MASS OF SPECIES
// --------------------------------------------------------------------------------------------

auto ChemicalState::setSpeciesAmounts(real const& value) -> void
{
    pimpl->setSpeciesAmounts(value);
}

auto ChemicalState::setSpeciesAmounts(ArrayXrConstRef const& n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::setSpeciesAmounts(ArrayXdConstRef const& n) -> void
{
    pimpl->setSpeciesAmounts(n);
}

auto ChemicalState::setSpeciesAmount(Index ispecies, real const& amount) -> void
{
    pimpl->setSpeciesAmount(ispecies, amount);
}

auto ChemicalState::setSpeciesAmount(StringOrIndex const& species, real amount, Chars unit) -> void
{
    pimpl->setSpeciesAmount(species, amount, unit);
}

auto ChemicalState::setSpeciesMass(StringOrIndex const& species, real mass, Chars unit) -> void
{
    pimpl->setSpeciesMass(species, mass, unit);
}

auto ChemicalState::set(StringOrIndex const& species, real value, Chars unit) -> void
{
    pimpl->set(species, value, unit);
}

auto ChemicalState::add(StringOrIndex const& species, real value, Chars unit) -> void
{
    pimpl->add(species, value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS FOR GETTING THE AMOUNT OR MASS OF SPECIES, ELEMENTS, AND CHARGE
// --------------------------------------------------------------------------------------------

auto ChemicalState::speciesAmounts() const -> ArrayXrConstRef
{
    return pimpl->n;
}

auto ChemicalState::speciesAmountsInPhase(StringOrIndex const& phase) const -> ArrayXrConstRef
{
    return pimpl->speciesAmountsInPhase(phase);
}

auto ChemicalState::speciesAmount(StringOrIndex const& species) const -> real
{
    return pimpl->speciesAmount(species);
}

auto ChemicalState::speciesMass(StringOrIndex const& species) const -> real
{
    return pimpl->speciesMass(species);
}

auto ChemicalState::componentAmounts() const -> ArrayXr
{
    return pimpl->componentAmounts();
}

auto ChemicalState::elementAmounts() const -> ArrayXr
{
    return pimpl->elementAmounts();
}

auto ChemicalState::charge() const -> real
{
    return pimpl->charge();
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE AMOUNTS OF SPECIES IN THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleSpeciesAmounts(real const& scalar) -> void
{
    pimpl->scaleSpeciesAmounts(scalar);
}

auto ChemicalState::scaleSpeciesAmounts(real const& scalar, Indices const& indices) -> void
{
    pimpl->scaleSpeciesAmounts(scalar, indices);
}

auto ChemicalState::scaleSpeciesAmountsInPhase(StringOrIndex const& phase, real const& scalar) -> void
{
    pimpl->scaleSpeciesAmountsInPhase(phase, scalar);
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE AMOUNT OF THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleAmount(real value, Chars unit) -> void
{
    pimpl->scaleAmount(value, unit);
}

auto ChemicalState::scalePhaseAmount(StringOrIndex const& phase, real value, Chars unit) -> void
{
    pimpl->scalePhaseAmount(phase, value, unit);
}

auto ChemicalState::scaleFluidAmount(real value, Chars unit) -> void
{
    pimpl->scaleFluidAmount(value, unit);
}

auto ChemicalState::scaleSolidAmount(real value, Chars unit) -> void
{
    pimpl->scaleSolidAmount(value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE MASS OF THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleMass(real value, Chars unit) -> void
{
    pimpl->scaleMass(value, unit);
}

auto ChemicalState::scalePhaseMass(StringOrIndex const& phase, real value, Chars unit) -> void
{
    pimpl->scalePhaseMass(phase, value, unit);
}

auto ChemicalState::scaleFluidMass(real value, Chars unit) -> void
{
    pimpl->scaleFluidMass(value, unit);
}

auto ChemicalState::scaleSolidMass(real value, Chars unit) -> void
{
    pimpl->scaleSolidMass(value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS TO SCALE THE VOLUME OF THE SYSTEM OR PART OF IT
// --------------------------------------------------------------------------------------------

auto ChemicalState::scaleVolume(real value, Chars unit) -> void
{
    pimpl->scaleVolume(value, unit);
}

auto ChemicalState::scalePhaseVolume(StringOrIndex const& phase, real value, Chars unit) -> void
{
    pimpl->scalePhaseVolume(phase, value, unit);
}

auto ChemicalState::scaleFluidVolume(real value, Chars unit) -> void
{
    pimpl->scaleFluidVolume(value, unit);
}

auto ChemicalState::scaleSolidVolume(real value, Chars unit) -> void
{
    pimpl->scaleSolidVolume(value, unit);
}

// --------------------------------------------------------------------------------------------
// METHODS FOR UPDATING CHEMICAL STATE AND ITS PROPERTIES
// --------------------------------------------------------------------------------------------

auto ChemicalState::update(real const& T, real const& P, ArrayXrConstRef const& n) -> void
{
    setTemperature(T);
    setPressure(P);
    setSpeciesAmounts(n);
    props().update(T, P, n);
}

auto ChemicalState::updateIdeal(real const& T, real const& P, ArrayXrConstRef const& n) -> void
{
    setTemperature(T);
    setPressure(P);
    setSpeciesAmounts(n);
    props().updateIdeal(T, P, n);
}

// --------------------------------------------------------------------------------------------
// MISCELLANEOUS METHODS
// --------------------------------------------------------------------------------------------

auto ChemicalState::system() const -> ChemicalSystem const&
{
    return pimpl->system;
}

auto ChemicalState::props() const -> ChemicalProps const&
{
    return pimpl->props;
}

auto ChemicalState::props() -> ChemicalProps&
{
    return pimpl->props;
}

auto ChemicalState::equilibrium() const -> Equilibrium const&
{
    return pimpl->equilibrium;
}

auto ChemicalState::equilibrium() -> Equilibrium&
{
    return pimpl->equilibrium;
}

auto ChemicalState::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ChemicalState::output(String const& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

//=================================================================================================
//
// ChemicalState::Equilibrium
//
//=================================================================================================

struct ChemicalState::Equilibrium::Impl
{
    /// The number of species in the chemical system.
    const Index Nn;

    /// The number of components in the equilibrium state.
    const Index Nb;

    /// The number of *q* control variables considered in the equilibrium calculation.
    Index Nq;

    /// The names of the *w* control variables considered in the equilibrium calculation.
    Strings wnames;

    /// The names of the *p* control variables considered in the equilibrium calculation.
    Strings pnames;

    /// The names of the *q* control variables considered in the equilibrium calculation.
    Strings qnames;

    /// The values of the input variables *w* used in the equilibrium calculation.
    ArrayXd w;

    /// The initial component amounts in the equilibrium calculation.
    ArrayXd c;

    /// The Optima::State object used for warm start Optima optimization calculations.
    Optima::State optstate;

    /// Construct a default ChemicalState::Equilibrium::Impl instance
    Impl(ChemicalSystem const& system)
    : Nn(system.species().size()), Nb(system.elements().size() + 1)
    {}
};

ChemicalState::Equilibrium::Equilibrium(ChemicalSystem const& system)
: pimpl(new Impl(system))
{}

ChemicalState::Equilibrium::Equilibrium(ChemicalState::Equilibrium const& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalState::Equilibrium::~Equilibrium()
{}

auto ChemicalState::Equilibrium::operator=(ChemicalState::Equilibrium other) -> Equilibrium&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalState::Equilibrium::reset() -> void
{
    pimpl->Nq = 0;
    pimpl->wnames = {};
    pimpl->pnames = {};
    pimpl->qnames = {};
    pimpl->w = {};
    pimpl->c = {};
    pimpl->optstate = {};
}

auto ChemicalState::Equilibrium::setNamesInputVariables(Strings const& wnames) -> void
{
    pimpl->wnames = wnames;
}

auto ChemicalState::Equilibrium::setNamesControlVariablesP(Strings const& pnames) -> void
{
    pimpl->pnames = pnames;
}

auto ChemicalState::Equilibrium::setNamesControlVariablesQ(Strings const& qnames) -> void
{
    pimpl->qnames = qnames;
    pimpl->Nq = qnames.size();
}

auto ChemicalState::Equilibrium::setInputVariables(ArrayXdConstRef const& w) -> void
{
    errorifnot(w.size() == pimpl->wnames.size(), "The number of input variables in the equilibrium calculation must be equal to the number of registered input variables in the specifications of the equilibrium problem. Have you constructed an EquilibriumSolver object with a ChemicalSystem object instead of an EquilibriumSpecs object (e.g., EquilibriumSolver(system) instead of EquilibriumSolver(specs))?");
    pimpl->w = w;
}

auto ChemicalState::Equilibrium::setControlVariablesP(ArrayXdConstRef const& p) -> void
{
    errorifnot(p.size() == pimpl->pnames.size(), "The number of input control variables p in the equilibrium calculation must be equal to the number of registered input control variables p in the specifications of the equilibrium problem. Have you constructed an EquilibriumSolver object with a ChemicalSystem object instead of an EquilibriumSpecs object (e.g., EquilibriumSolver(system) instead of EquilibriumSolver(specs))?");
    pimpl->optstate.p = p;
}

auto ChemicalState::Equilibrium::setControlVariablesQ(ArrayXdConstRef const& q) -> void
{
    errorifnot(q.size() == pimpl->qnames.size(), "The number of input control variables q in the equilibrium calculation must be equal to the number of registered input control variables q in the specifications of the equilibrium problem. Have you constructed an EquilibriumSolver object with a ChemicalSystem object instead of an EquilibriumSpecs object (e.g., EquilibriumSolver(system) instead of EquilibriumSolver(specs))?");
    if(pimpl->Nq > 0)
        pimpl->optstate.x.tail(pimpl->Nq) = q;
}

auto ChemicalState::Equilibrium::setInitialComponentAmounts(ArrayXdConstRef const& c) -> void
{
    pimpl->c = c;
}

auto ChemicalState::Equilibrium::setOptimaState(Optima::State const& state) -> void
{
    pimpl->optstate = state;
}

auto ChemicalState::Equilibrium::empty() const -> bool
{
    return pimpl->optstate.x.size() == 0; // this means optstate has not been set yet
}

auto ChemicalState::Equilibrium::numPrimarySpecies() const -> Index
{
    return pimpl->optstate.jb.size();
}

auto ChemicalState::Equilibrium::numSecondarySpecies() const -> Index
{
    return pimpl->optstate.jn.size();
}

auto ChemicalState::Equilibrium::indicesPrimarySpecies() const -> ArrayXlConstRef
{
    return pimpl->optstate.jb;
}

auto ChemicalState::Equilibrium::indicesSecondarySpecies() const -> ArrayXlConstRef
{
    return pimpl->optstate.jn;
}

auto ChemicalState::Equilibrium::elementChemicalPotentials() const -> ArrayXdConstRef
{
    if(pimpl->optstate.ye.size())
        return pimpl->optstate.ye.head(pimpl->Nb);
    else return pimpl->optstate.ye;
}

auto ChemicalState::Equilibrium::speciesStabilities() const -> ArrayXdConstRef
{
    if(pimpl->optstate.s.size())
        return pimpl->optstate.s.head(pimpl->Nn);
    else return pimpl->optstate.s;
}

auto ChemicalState::Equilibrium::explicitTitrantAmount(String const& name) const -> real
{
    const auto idx = index(pimpl->pnames, "[" + name + "]");
    errorif(idx >= pimpl->pnames.size(), "There is no explicit titrant with name `", name, "` in this ChemicalState object.");
    return -pimpl->optstate.p[idx]; // note negative sign due to convention for p when used for titrant amounts so that the conservation matrix has positive element coefficients for the explicit titrants, which showed to work better algorithmically
}

auto ChemicalState::Equilibrium::implicitTitrantAmount(String const& name) const -> real
{
    const auto idx = index(pimpl->qnames, "[" + name + "]");
    errorif(idx >= pimpl->qnames.size(), "There is no implicit titrant with name `", name, "` in this ChemicalState object.");
    return -pimpl->optstate.x[pimpl->Nn + idx]; // note negative sign due to convention for q so that the conservation matrix has positive element coefficients for the implicit titrants, which showed to work better algorithmically
}

auto ChemicalState::Equilibrium::titrantAmount(String const& name) const -> real
{
    const auto pidx = index(pimpl->pnames, "[" + name + "]");
    if(pidx < pimpl->pnames.size())
        return -pimpl->optstate.p[pidx]; // note negative sign due to convention for p when used for titrant amounts so that the conservation matrix has positive element coefficients for the explicit titrants, which showed to work better algorithmically

    const auto qidx = index(pimpl->qnames, "[" + name + "]");
    if(qidx < pimpl->qnames.size())
        return -pimpl->optstate.x[pimpl->Nn + qidx]; // note negative sign due to convention for q so that the conservation matrix has positive element coefficients for the implicit titrants, which showed to work better algorithmically

    errorif(true, "There is no explicit nor implicit titrant with name `", name, "` in this ChemicalState object.");
}

auto ChemicalState::Equilibrium::namesInputVariables() const -> Strings const&
{
    return pimpl->wnames;
}

auto ChemicalState::Equilibrium::namesControlVariablesP() const -> Strings const&
{
    return pimpl->pnames;
}

auto ChemicalState::Equilibrium::namesControlVariablesQ() const -> Strings const&
{
    return pimpl->qnames;
}

auto ChemicalState::Equilibrium::inputVariables() const -> ArrayXdConstRef
{
    return pimpl->w;
}

auto ChemicalState::Equilibrium::controlVariablesP() const -> ArrayXdConstRef
{
    return pimpl->optstate.p;
}

auto ChemicalState::Equilibrium::controlVariablesQ() const -> ArrayXdConstRef
{
    return pimpl->optstate.x.tail(pimpl->Nq);
}

auto ChemicalState::Equilibrium::initialComponentAmounts() const -> ArrayXdConstRef
{
    return pimpl->c;
}

auto ChemicalState::Equilibrium::w() const -> ArrayXdConstRef
{
    return inputVariables();
}

auto ChemicalState::Equilibrium::p() const -> ArrayXdConstRef
{
    return controlVariablesP();
}

auto ChemicalState::Equilibrium::q() const -> ArrayXdConstRef
{
    return controlVariablesQ();
}

auto ChemicalState::Equilibrium::c() const -> ArrayXdConstRef
{
    return initialComponentAmounts();
}

auto ChemicalState::Equilibrium::optimaState() const -> Optima::State const&
{
    return pimpl->optstate;
}

auto operator<<(std::ostream& out, ChemicalState const& state) -> std::ostream&
{
    auto const& n = state.speciesAmounts();
    auto const& b = state.elementAmounts();
    auto const& phases = state.system().phases();
    auto const& species = state.system().species();
    auto const& elements = state.system().elements();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(state.temperature()), "K" });
    table.add_row({ "Pressure", strfix(state.pressure()*1e-5), "bar" });
    table.add_row({ "Charge:", strsci(state.charge()), "mol" });

    table.add_row({ "Element Amount:", "", "" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), strsci(b[i]), "mol" });
    table.add_row({ "Species Amount:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(n[i]), "mol" });

    auto i = 0;
    for(auto& row : table)
    {
        if(i >= 2)  // apply from the third row
            table[i]
                .format()
                .border_top("")
                .column_separator("")
                .corner_top_left("")
                .corner_top_right("");
        i += 1;
    }

    table.row(0).format().font_style({FontStyle::bold});  // Bold face for header
    table.column(1).format().font_align(FontAlign::right); // Value column with right alignment
    table.column(2).format().font_align(FontAlign::right); // Unit column with right alignment

    auto old_locale = std::locale::global(std::locale("C")); // This locale logic is needed to avoid UnicodeDecodeError: 'utf-8' codec can't decode byte 0xa0 in position ...
    out << table;
    std::locale::global(old_locale);

    return out;
}

} // namespace Reaktoro
