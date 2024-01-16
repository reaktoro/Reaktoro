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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;
class ChemicalSystem;
class EquilibriumConditions;
class EquilibriumRestrictions;
class EquilibriumSpecs;
class KineticsSensitivity;
struct SmartKineticsOptions;
struct SmartKineticsResult;

/// Used for chemical kinetics calculations.
class SmartKineticsSolver
{
public:
    /// Construct an SmartKineticsSolver object with given chemical system.
    explicit SmartKineticsSolver(ChemicalSystem const& system);

    /// Construct an SmartKineticsSolver object with given chemical equilibrium specifications to be attained during chemical kinetics.
    explicit SmartKineticsSolver(EquilibriumSpecs const& specs);

    /// Construct a copy of an SmartKineticsSolver object.
    SmartKineticsSolver(SmartKineticsSolver const& other);

    /// Destroy this SmartKineticsSolver object.
    ~SmartKineticsSolver();

    /// Assign a copy of an SmartKineticsSolver object to this.
    auto operator=(SmartKineticsSolver other) -> SmartKineticsSolver&;

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS
    //
    //=================================================================================================================

    /// React a chemical state for a given time interval.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param dt The time step in the kinetics calculation (in s).
    auto solve(ChemicalState& state, real const& dt) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, real const& dt, EquilibriumRestrictions const& restrictions) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given constraint conditions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param conditions The specified constraint conditions to be attained during chemical kinetics
    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given constraint conditions and reactivity restrictions.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param conditions The specified constraint conditions to be attained during chemical kinetics
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartKineticsResult;

    //=================================================================================================================
    //
    // CHEMICAL KINETICS METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    /// React a chemical state for a given time interval and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the reacted state with respect to given input conditions
    /// @param dt The time step in the kinetics calculation (in s).
    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the reacted state with respect to given input conditions
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumRestrictions const& restrictions) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given constraint conditions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the reacted state with respect to given input conditions
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param conditions The specified constraint conditions to be attained during chemical kinetics
    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions) -> SmartKineticsResult;

    /// React a chemical state for a given time interval respecting given constraint conditions and reactivity restrictions and compute sensitivity derivatives.
    /// @param[in,out] state The initial guess for the calculation (in) and the computed reacted state (out)
    /// @param[out] sensitivity The sensitivity derivatives of the reacted state with respect to given input conditions
    /// @param dt The time step in the kinetics calculation (in s).
    /// @param conditions The specified constraint conditions to be attained during chemical kinetics
    /// @param restrictions The reactivity restrictions on the amounts of selected species
    auto solve(ChemicalState& state, KineticsSensitivity& sensitivity, real const& dt, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartKineticsResult;

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the kinetics solver.
    auto setOptions(SmartKineticsOptions const& options) -> void;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
