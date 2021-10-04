// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
class ChemicalSystem;
class PhaseList;
class ElementList;
class SpeciesList;

namespace detail {

/// Return the molar masses of the species.
auto molarMasses(const SpeciesList& species) -> ArrayXd;

/// Compute the amount of a species given a value in mass or amount unit.
/// @param system The chemical system in which the species is.
/// @param ispecies The index of the species in the chemical system.
/// @param value The quantity value of the species.
/// @param unit The quantity unit of the species.
auto computeSpeciesAmount(const ChemicalSystem& system, Index ispecies, real value, const String& unit) -> real;

/// Resolve the index of a elementn in a chemical system with given either aa element symbol or its index itself.
auto resolveElementIndex(const ChemicalSystem& system, StringOrIndex element) -> Index;

/// Resolve the index of a species in a chemical system with given either a species name or its index itself.
auto resolveSpeciesIndex(const ChemicalSystem& system, StringOrIndex species) -> Index;

/// Resolve the index of a phase in a chemical system with given either a phase name or its index itself.
auto resolvePhaseIndex(const ChemicalSystem& system, StringOrIndex phase) -> Index;

/// Assemble the formula matrix of the list of `species` with respect to given `elements`.
auto assembleFormulaMatrix(const SpeciesList& species, const ElementList& elements) -> MatrixXd;

// Extract names of the species from the species' list
auto extractNames(const SpeciesList& list) -> Strings;

// Extract names of the elements from the elements' list
auto extractNames(const ElementList& list) -> Strings;

// Extract names of the phases from the phases' list
auto extractNames(const PhaseList& list) -> Strings;

} // namespace detail
} // namespace Reaktoro
