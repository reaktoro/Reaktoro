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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used to describe the state of an aqueous mixture.
/// @see AqueousMixture
struct AqueousMixtureState
{
    /// The temperature of the mixture (in K).
    real T;

    /// The pressure of the mixture (in Pa).
    real P;

    /// The density of water (in kg/m3)
    real rho;

    /// The relative dielectric constant of water (no units)
    real epsilon;

    /// The effective ionic strength of the aqueous mixture (in mol/kg)
    real Ie;

    /// The stoichiometric ionic strength of the aqueous mixture (in mol/kg)
    real Is;

    /// The molalities of the aqueous species (in mol/kg)
    ArrayXr m;

    /// The stoichiometric molalities of the ionic species (in mol/kg)
    ArrayXr ms;
};

/// A type used to describe an aqueous mixture.
/// The AqueousMixture class is defined as a collection of Species objects,
/// representing, therefore, a mixture of aqueous species. Its main purpose is to
/// provide the necessary operations in the calculation of activities of aqueous
/// species. It implements methods for the calculation of mole fractions, molalities,
/// stoichiometric molalities, and effective and stoichiometric ionic strengths.
/// In addition, it provides methods that retrieves information about the ionic, neutral
/// and complex species.
/// @ingroup GeochemistryExtension
class AqueousMixture
{
public:
    /// Construct a default AqueousMixture instance.
    AqueousMixture();

    /// Construct an AqueousMixture instance with given species.
    explicit AqueousMixture(SpeciesList const& species);

    /// Return a deep copy of this AqueousMixture object.
    auto clone() const -> AqueousMixture;

    /// Return a copy of this AqueousMixture object with replaced function for water density calculation.
    auto withWaterDensityFn(Fn<real(real,real)> rho) const -> AqueousMixture;

    /// Return a copy of this AqueousMixture object with replaced function for water dielectric constant calculation.
    auto withWaterDielectricConstantFn(Fn<real(real,real)> epsilon) const -> AqueousMixture;

    /// Return the aqueous species in the mixture with given index.
    auto species(Index idx) const -> Species const&;

    /// Return the aqueous species in the mixture.
    auto species() const -> SpeciesList const&;

    /// Return the neutral aqueous solutes in the mixture.
    auto neutral() const -> SpeciesList const&;

    /// Return the charged aqueous solutes in the mixture.
    auto charged() const -> SpeciesList const&;

    /// Return the cation solutes in the mixture.
    auto cations() const -> SpeciesList const&;

    /// Return the anion solutes in the mixture.
    auto anions() const -> SpeciesList const&;

    /// Return the aqueous solvent species in the mixture.
    auto water() const -> Species const&;

    /// Return the indices of the neutral aqueous solutes in the mixture.
    auto indicesNeutral() const -> Indices const&;

    /// Return the indices of the charged aqueous solutes in the mixture.
    auto indicesCharged() const -> Indices const&;

    /// Return the indices of the cations in the mixture.
    auto indicesCations() const -> Indices const&;

    /// Return the indices of the anions in the mixture.
    auto indicesAnions() const -> Indices const&;

    /// Return the index of the solvent species in the mixture.
    auto indexWater() const -> Index;

    /// Return the electric charges of the aqueous species in the mixture.
    auto charges() const -> ArrayXdConstRef;

    /// Return the dissociation matrix of the neutral species into charged species.
    /// The dissociation matrix of the aqueous mixture is defined so that its
    /// entry (i, j) corresponds to the stoichiometric coefficient of the j-th
    /// charged species in the i-th neutral species. It is used to compute
    /// stoichiometric molalities of the charged species as well as the
    /// stoichiometric ionic strength of the mixture.
    auto dissociationMatrix() const -> MatrixXdConstRef;

    /// Calculate the state of the aqueous mixture.
    /// @param T The temperature (in K)
    /// @param P The pressure (in Pa)
    /// @param x The mole fractions of the species in the mixture
    auto state(real T, real P, ArrayXrConstRef x) const -> AqueousMixtureState;

    /// Set the default function for water density calculation when creating AqueousMixture objects.
    static auto setDefaultWaterDensityFn(Fn<real(real,real)> rho) -> void;

    /// Set the default function for water dielectric constant calculation when creating AqueousMixture objects.
    static auto setDefaultWaterDielectricConstantFn(Fn<real(real,real)> epsilon) -> void;

    /// Reset the default function for water density calculation when creating AqueousMixture objects.
    /// The default function is the one that returns the density of water at 25°C and 1 atm:
    /// ρ = 997.04703901770279 kg/m3 (via Wagner-Pruss equation)
    static auto resetDefaultWaterDensityFn() -> void;

    /// Reset the default function for water dielectric constant calculation when creating AqueousMixture objects.
    /// The default function is the one that returns the density of water at 25°C and 1 atm:
    /// ϵ = 78.245144808202397 (via Jonhson-Norton equation)
    static auto resetDefaultWaterDielectricConstantFn() -> void;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

} // namespace Reaktoro
