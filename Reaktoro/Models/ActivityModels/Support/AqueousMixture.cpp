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

#include "AqueousMixture.hpp"

// C++ includes
#include <algorithm>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Singletons/DissociationReactions.hpp>
#include <Reaktoro/Water/WaterElectroProps.hpp>
#include <Reaktoro/Water/WaterElectroPropsJohnsonNorton.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>
#include <Reaktoro/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace detail {

auto defaultWaterDensityFn() -> Fn<real(real,real)>
{
    const real rho0 = 997.04703901770279; // density value in kg/m3 previously obtained by executing waterLiquidDensityWagnerPruss(T, P) with T = 298.15 K (25 °C) and P = 100,000 Pa (1 bar)
    return [=](real T, real P) { return rho0; };
}

auto defaultWaterDielectricConstantFn() -> Fn<real(real,real)>
{
    const real epsilon0 = 78.245144808202397; // dielectric constant (no unit) previously obtained by executing waterElectroPropsJohnsonNorton(T, P, waterThermoPropsHGK(T, P, StateOfMatter::Liquid)) with T = 298.15 K (25 °C) and P = 100,000 Pa (1 bar)
    return [=](real T, real P) { return epsilon0; };
}

auto default_water_density_fn = defaultWaterDensityFn();
auto default_water_dielectric_constant_fn = defaultWaterDielectricConstantFn();

} // namespace detail

struct AqueousMixture::Impl
{
    /// The aqueous species in the mixture.
    SpeciesList species;

    /// The neutral aqueous species in the mixture.
    SpeciesList neutral;

    /// The charged aqueous species in the mixture.
    SpeciesList charged;

    /// The cation species in the mixture.
    SpeciesList cations;

    /// The anion species in the mixture.
    SpeciesList anions;

    /// The aqueous solvent species in the mixture.
    Species water;

    /// The indices of the neutral aqueous species in the mixture.
    Indices idx_neutral_species;

    /// The indices of the charged aqueous species in the mixture.
    Indices idx_charged_species;

    /// The indices of the cation species in the mixture.
    Indices idx_cations;

    /// The indices of the anion species in the mixture.
    Indices idx_anions;

    /// The index of the aqueous solvent species in the mixture.
    Index idx_water;

    /// The electric charges of the aqueous species in the mixture.
    ArrayXd z;

    /// The matrix that represents the dissociation of the aqueous complexes into ions.
    MatrixXd dissociation_matrix;

    /// The density function for water.
    Fn<real(real,real)> rho;

    /// The dielectric constant function for water.
    Fn<real(real,real)> epsilon;

    /// Construct a default AqueousMixture::Impl instance.
    Impl()
    {}

    /// Construct an AqueousMixture::Impl instance with given species.
    Impl(SpeciesList const& species)
    : species(species)
    {
        // Initialize the index related data
        initializeIndices();

        /// Initialize the array of electric charges of the species
        initializeCharges();

        // Initialize the dissociation matrix of the neutral species w.r.t. the charged species
        initializeDissociationMatrix();

        // Initialize the density function for water
        rho = detail::default_water_density_fn;

        // Initialize the dielectric constant function for water
        epsilon = detail::default_water_dielectric_constant_fn;
    }

    /// Initialize the index related data of the species.
    auto initializeIndices() -> void
    {
        // Initialize the water species and its index
        idx_water = species.indexWithFormula("H2O");
        water = species[idx_water];

        // Initialize the indices of the charged and neutral species
        for(auto i = 0; i < species.size(); ++i)
        {
            if(i == idx_water)
                continue;
            if(species[i].charge() == 0.0)
            {
                idx_neutral_species.push_back(i);
                neutral.push_back(species[i]);
            }
            else
            {
                idx_charged_species.push_back(i);
                charged.push_back(species[i]);
                if(species[i].charge() > 0.0)
                {
                    idx_cations.push_back(i);
                    cations.push_back(species[i]);
                }
                else
                {
                    idx_anions.push_back(i);
                    anions.push_back(species[i]);
                }
            }
        }
    }

    /// Initialize the array of electric charges of the species
    auto initializeCharges() -> void
    {
        const auto charges = vectorize(species, RKT_LAMBDA(x, x.charge()));
        z = ArrayXd::Map(charges.data(), charges.size());
    }

    /// Initialize the dissociation matrix of the neutral species w.r.t. the charged species.
    auto initializeDissociationMatrix() -> void
    {
        // Return the stoichiometry of the i-th charged species in the j-th neutral species
        auto stoichiometry = [&](Index i, Index j)
        {
            const auto ineutral = idx_neutral_species[i];
            const auto icharged = idx_charged_species[j];
            const auto neutral = species[ineutral].formula();
            const auto charged = species[icharged].formula();
            return DissociationReactions::coefficient(neutral, charged);
        };

        // Assemble the dissociation matrix of the neutral species with respect to the charged species
        const auto num_charged_species = idx_charged_species.size();
        const auto num_neutral_species = idx_neutral_species.size();
        dissociation_matrix.resize(num_neutral_species, num_charged_species);
        for(auto i = 0; i < num_neutral_species; ++i)
            for(auto j = 0; j < num_charged_species; ++j)
                dissociation_matrix(i, j) = stoichiometry(i, j);
    }

    /// Return the molalities of the aqueous species with given mole fractions.
    auto molalities(ArrayXrConstRef x) const -> ArrayXr
    {
        const auto xw = x[idx_water];
        const auto Mw = water.molarMass();
        if(xw == 0.0)
            return ArrayXr::Zero(x.size());
        return x/(Mw * xw);
    }

    /// Return the stoichiometric molalities of the charged species with given molalities.
    auto stoichiometricMolalities(ArrayXrConstRef m) const -> ArrayXr
    {
        // The molalities of the charged species
        const auto mc = m(idx_charged_species).matrix(); // convert from array to matrix expression

        // The molalities of the neutral species
        const auto mn = m(idx_neutral_species).matrix(); // convert from array to matrix expression

        // The stoichiometric molalities of the charged species
        const auto ms = mc + dissociation_matrix.transpose() * mn;

        return ms;
    }

    /// Return the effective ionic strength of the aqueous mixture with given molalities.
    auto effectiveIonicStrength(ArrayXrConstRef m) const -> real
    {
        return 0.5 * (z * z * m).sum();
    }

    /// Return the stoichiometric ionic strength of the aqueous mixture with given stoichiometric molalities of the charged species.
    auto stoichiometricIonicStrength(ArrayXrConstRef ms) const -> real
    {
        const auto zc = z(idx_charged_species);
        return 0.5 * (zc * zc * ms).sum();
    }

    /// Return the state of the aqueous mixture.
    auto state(real T, real P, ArrayXrConstRef x) const -> AqueousMixtureState
    {
        AqueousMixtureState state;
        state.T = T;
        state.P = P;
        state.rho = rho(T, P);
        state.epsilon = epsilon(T, P);
        state.m  = molalities(x);
        state.ms = stoichiometricMolalities(state.m);
        state.Ie = effectiveIonicStrength(state.m);
        state.Is = stoichiometricIonicStrength(state.ms);
        return state;
    }
};

AqueousMixture::AqueousMixture()
: pimpl(new Impl())
{}

AqueousMixture::AqueousMixture(SpeciesList const& species)
: pimpl(new Impl(species))
{}

auto AqueousMixture::clone() const -> AqueousMixture
{
    AqueousMixture copy;
    *copy.pimpl = *pimpl;
    return copy;
}

auto AqueousMixture::withWaterDensityFn(Fn<real(real,real)> rho) const -> AqueousMixture
{
    AqueousMixture copy = clone();
    copy.pimpl->rho = std::move(rho);
    return copy;
}

auto AqueousMixture::withWaterDielectricConstantFn(Fn<real(real,real)> epsilon) const -> AqueousMixture
{
    AqueousMixture copy = clone();
    copy.pimpl->epsilon = std::move(epsilon);
    return copy;
}

auto AqueousMixture::species(Index idx) const -> Species const&
{
    return pimpl->species[idx];
}

auto AqueousMixture::species() const -> SpeciesList const&
{
    return pimpl->species;
}

auto AqueousMixture::neutral() const -> SpeciesList const&
{
    return pimpl->neutral;
}

auto AqueousMixture::charged() const -> SpeciesList const&
{
    return pimpl->charged;
}

auto AqueousMixture::cations() const -> SpeciesList const&
{
    return pimpl->cations;
}

auto AqueousMixture::anions() const -> SpeciesList const&
{
    return pimpl->anions;
}

auto AqueousMixture::water() const -> Species const&
{
    return pimpl->water;
}

auto AqueousMixture::indicesNeutral() const -> Indices const&
{
    return pimpl->idx_neutral_species;
}

auto AqueousMixture::indicesCharged() const -> Indices const&
{
    return pimpl->idx_charged_species;
}

auto AqueousMixture::indicesCations() const -> Indices const&
{
    return pimpl->idx_cations;
}

auto AqueousMixture::indicesAnions() const -> Indices const&
{
    return pimpl->idx_anions;
}

auto AqueousMixture::indexWater() const -> Index
{
    return pimpl->idx_water;
}

auto AqueousMixture::dissociationMatrix() const -> MatrixXdConstRef
{
    return pimpl->dissociation_matrix;
}

auto AqueousMixture::charges() const -> ArrayXdConstRef
{
    return pimpl->z;
}

auto AqueousMixture::state(real T, real P, ArrayXrConstRef x) const -> AqueousMixtureState
{
    return pimpl->state(T, P, x);
}

auto AqueousMixture::setDefaultWaterDensityFn(Fn<real(real,real)> rho) -> void
{
    detail::default_water_density_fn = std::move(rho);
}

auto AqueousMixture::setDefaultWaterDielectricConstantFn(Fn<real(real,real)> epsilon) -> void
{
    detail::default_water_dielectric_constant_fn = std::move(epsilon);
}

auto AqueousMixture::resetDefaultWaterDensityFn() -> void
{
    detail::default_water_density_fn = detail::defaultWaterDensityFn();
}

auto AqueousMixture::resetDefaultWaterDielectricConstantFn() -> void
{
    detail::default_water_dielectric_constant_fn = detail::defaultWaterDielectricConstantFn();
}

} // namespace Reaktoro
