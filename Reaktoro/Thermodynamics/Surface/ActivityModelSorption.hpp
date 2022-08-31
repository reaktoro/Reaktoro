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
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Thermodynamics/Surface/Surface.hpp>

namespace Reaktoro {

// The parameters in the Diffuse Double Layer (DDL) model for the diffusive double layer solution between complexation
/// surface and aqueous solution.
/// @ingroup ActivityModels
struct ActivityModelDDLParams
{
    /// The default value of the enrichment factor in the Donnan pore space.
    real enr = 1e+0; // TODO: influence the results quite a lot, figure out why
    // (e.g., enr = 1e-2, puts all the amounts to the DDL species)

    /// The diffusive double layer thickness.
    real thickness = 1e-8;

    // Output option
    bool output = false;
};

/// The parameters for the surface complexation activity model (considering every site as a separate phase).
/// @ingroup ActivityModels
struct ActivityModelSorptionParams
{
    /// The complexation surface.
    Surface surface;

    String site_tag;

    // The double layer.
    ActivityModelDDLParams ddl;

    // Output option
    bool output = false;
};

///// The parameters in the Diffuse Double Layer (DDL) model for the diffusive double layer solution between complexation
///// surface and aqueous solution.
///// @ingroup ActivityModels
//struct ActivityModelDDLDonnanParams
//{
//    ActivityModelDDLParams ddl;
//
//    /// The Debye length of the diffusive double layer.
//    real debye_lengths = 1.0;
//
//    /// The amount of water contained in diffusive layer can be limited.
//    real limit = 0.8; // the fraction of the total water (pore space plus diffuse double layer water) that can be in the diffuse double layer
//
//    /// The viscosity.
//    real viscosity = 1.0;
//};

/// Return the SorptionActivityModel object assuming no electrostatic effects are generated by the surface site.
auto ActivityModelSorptionNoDDL(ActivityModelSorptionParams params) -> ActivityModelGenerator;

///// Return the SorptionActivityModel object updating the surface potential of the surface site (must be combined
///// with with activity model applied on the aqueous solution that is being influenced by these electrostatic effects).
//auto ActivityModelSorptionWithDDL(ActivityModelSorptionParams params) -> ActivityModelGenerator;
//
///// Return the SorptionActivityModel object updating the surface potential of the surface site and applying
///// the electrostatic (coloumbic) corrections on the surface site species based on the Dzombak and Morel (1990) book.
//auto ActivityModelSorptionWithEDL(ActivityModelSorptionParams params) -> ActivityModelGenerator;
//
///// Return the activity model for the Aqueous phase, which include electrostatic correction of the ions.
///// @ingroup ActivityModels
//auto ActivityModelEDL(ActivityModelDDLParams params) -> ActivityModelGenerator;

//=====================================================================================================================
/// @page PageActivityModelSorptionDDL Diffuse Double Layer (DDL) activity model
///
/// The Diffuse Double Layer (DDL) model (from Dzombak & Morel, 1990), the Gouy-Chapman Diffuse Layer.
/// Apparent equilibrium constant @eq{K_{app}} (the value of which is dependent on the surface charge) is calculated
/// with a *coulombic correction* used by Dzombak and Morel (1990), i.e.,
/// @eqc{K_{app} = K_{int}\exp{\dfrac{{\Delta}ZF\Psi}{RT}}, }
/// where @eq{\Delta Z} is the change in the charge of the surface species due to
/// sorption reactions, @eq{\Psi} is the surface potential in volts, F is the Faraday
/// constant (96,485 C/mol), R is the molar gas constant (8.314 J/(mol*K))
/// and T is the absolute temperature in K.
///
/// The surface potential @eq{\Psi} (volt) is defined via the charge density @eq{\sigma} (C/m2) via the relation
/// @eqc{\sigma = 0.1174 I^{1⁄2} \sinh(\dfrac{F\Psi}{2RT}),}
/// where @eq{I} is the ionic strength.
///
/// The charge density of the surface is denoted by @eqc{\sigma} (C/m2) and defined as
/// @eqc{\sigma = \dfrac{FZ}{A_s m},}
/// where
/// @eq{F} is the Faraday constant (96,485 C/mol),
/// @eq{Z = \sum_{j}x_{j}z_{j}}}} is the surface charge (eq), defined as sum of the species molalities multiplied
/// by the the corresponding charge,
/// @eq{A_s} is specific surface area (m2/kg),
/// @eq{m} is the solid concentration in (kg).
///
//=====================================================================================================================


} // namespace Reaktoro