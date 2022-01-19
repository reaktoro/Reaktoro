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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ReactionThermoProps.hpp>

namespace Reaktoro {

/// The parameters in a thermodynamic model for a formation reaction based on constant @eq{\lg K(T)}.
struct ReactionThermoModelParamsConstLgK
{
    /// The equilibrium constant @eq{\lg K_{\mathrm{r}}} (log base 10) of the reaction at @eq{T_{\mathrm{r}}} and @eq{P_{\mathrm{r}}}.
    Param lgKr;

    /// The reference pressure @eq{P_{\mathrm{r}}} (in Pa).
    real Pr = 100'000;
};

/// Return a function that calculates thermodynamic properties of a reaction using a constant model for @eq{\lg K(T)}.
///
/// In this model, the equilibrium constant of the reaction remains constant:
///
/// @eqc{\lg K(T)=\lg K_{\mathrm{r}},}
///
/// where @eq{\lg K_{\mathrm{r}}} is a given equilibrium constant at a
/// reference temperature and pressure.
///
/// The standard Gibbs energy of the reaction is computed using:
///
/// @eqc{\Delta G^{\circ}(T,P)=-RT\ln K_{\mathrm{r}}+\Delta V^{\circ}(P-P_{\mathrm{r}})}
///
/// while the standard enthalpy of the reaction is:
///
/// @eqc{\Delta H^{\circ}(T,P)=\Delta V^{\circ}(P-P_{\mathrm{r}}).}
///
/// Note that a pressure correction is introduced above, where @eq{\Delta V^{\circ}}
/// is the change of standard molar volume of the reaction at @eq{T} and @eq{P},
/// with @eq{P_{\mathrm{r}}} denoting a given reference pressure.
auto ReactionThermoModelConstLgK(const ReactionThermoModelParamsConstLgK& params) -> ReactionThermoModel;

} // namespace Reaktoro
