// Reaktoro is a unified framework for modeling chemically reactive phases.
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

#include "ChemicalProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

ChemicalProps::ChemicalProps()
{}

ChemicalProps::ChemicalProps(ChemicalSystem const& system)
: msystem(system)
{
    const auto N = system.species().size();
    const auto K = system.phases().size();

    Ts   = ArrayXr::Zero(K);
    Ps   = ArrayXr::Zero(K);
    n    = ArrayXr::Zero(N);
    nsum = ArrayXr::Zero(K);
    msum = ArrayXr::Zero(K);
    x    = ArrayXr::Zero(N);
    G0   = ArrayXr::Zero(N);
    H0   = ArrayXr::Zero(N);
    V0   = ArrayXr::Zero(N);
    VT0  = ArrayXr::Zero(N);
    VP0  = ArrayXr::Zero(N);
    Cp0  = ArrayXr::Zero(N);
    Vx   = ArrayXr::Zero(K);
    VxT  = ArrayXr::Zero(K);
    VxP  = ArrayXr::Zero(K);
    Vxi  = ArrayXr::Zero(N);
    Gx   = ArrayXr::Zero(K);
    Hx   = ArrayXr::Zero(K);
    Cpx  = ArrayXr::Zero(K);
    ln_g = ArrayXr::Zero(N);
    ln_a = ArrayXr::Zero(N);
    u    = ArrayXr::Zero(N);
    som.resize(K);
}

ChemicalProps::ChemicalProps(ChemicalState const& state)
: ChemicalProps(state.system())
{
    update(state);
}

auto ChemicalProps::update(ChemicalState const& state) -> void
{
    auto const& T = state.temperature();
    auto const& P = state.pressure();
    auto const& n = state.speciesAmounts();
    update(T, P, n);
}

auto ChemicalProps::update(real const& T0, real const& P0, ArrayXrConstRef n0) -> void
{
    mstateid += 1;

    assert(T0 >= 0.0);
    assert(P0 >= 0.0);
    assert(n0.size() == n.size() && (n0 >= 0.0).all());

    T = T0;
    P = P0;

    auto offset = 0;
    for(auto const& [i, phase] : enumerate(msystem.phases()))
    {
        const auto size = phase.species().size();
        const auto np = n0.segment(offset, size);
        phasePropsRef(i).update(T, P, np, m_extra);
        offset += size;
    }
}

auto ChemicalProps::update(ArrayXrConstRef data) -> void
{
    mstateid += 1;
    ArraySerialization::deserialize(data, T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::update(ArrayXdConstRef data) -> void
{
    mstateid += 1;
    ArraySerialization::deserialize(data, T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::updateIdeal(ChemicalState const& state) -> void
{
    auto const& T = state.temperature();
    auto const& P = state.pressure();
    auto const& n = state.speciesAmounts();
    updateIdeal(T, P, n);
}

auto ChemicalProps::updateIdeal(real const& T0, real const& P0, ArrayXrConstRef n0) -> void
{
    mstateid += 1;

    assert(T0 >= 0.0);
    assert(P0 >= 0.0);
    assert(n0.size() == n.size() && (n0 >= 0.0).all());

    T = T0;
    P = P0;

    auto offset = 0;
    for(auto const& [i, phase] : enumerate(msystem.phases()))
    {
        const auto size = phase.species().size();
        const auto np = n0.segment(offset, size);
        phasePropsRef(i).updateIdeal(T, P, np, m_extra);
        offset += size;
    }
}

auto ChemicalProps::serialize(ArrayStream<real>& stream) const -> void
{
    stream.from(T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::serialize(ArrayStream<double>& stream) const -> void
{
    stream.from(T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::deserialize(const ArrayStream<real>& stream) -> void
{
    mstateid += 1;
    stream.to(T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::deserialize(const ArrayStream<double>& stream) -> void
{
    mstateid += 1;
    stream.to(T, P, n, Ts, Ps, nsum, msum, x, G0, H0, V0, VT0, VP0, Cp0, Vx, VxT, VxP, Vxi, Gx, Hx, Cpx, ln_g, ln_a, u);
}

auto ChemicalProps::stateid() const -> Index
{
    return mstateid;
}

auto ChemicalProps::system() const -> ChemicalSystem const&
{
    return msystem;
}

auto ChemicalProps::phaseProps(StringOrIndex phaseid) const -> ChemicalPropsPhaseConstRef
{
    return const_cast<ChemicalProps&>(*this).phasePropsRef(phaseid);
}

auto ChemicalProps::phasePropsRef(StringOrIndex phaseid) -> ChemicalPropsPhaseRef
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phaseid);
    const auto phase = msystem.phase(iphase);
    const auto begin = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto size = phase.species().size();

    return ChemicalPropsPhaseRef(phase, {
        Ts[iphase],
        Ps[iphase],
        n.segment(begin, size),
        nsum[iphase],
        msum[iphase],
        x.segment(begin, size),
        G0.segment(begin, size),
        H0.segment(begin, size),
        V0.segment(begin, size),
        VT0.segment(begin, size),
        VP0.segment(begin, size),
        Cp0.segment(begin, size),
        Vx[iphase],
        VxT[iphase],
        VxP[iphase],
        Vxi.segment(begin, size),
        Gx[iphase],
        Hx[iphase],
        Cpx[iphase],
        ln_g.segment(begin, size),
        ln_a.segment(begin, size),
        u.segment(begin, size),
        som[iphase],
    });
}

auto ChemicalProps::extra() const -> const Map<String, Any>&
{
    return m_extra;
}

auto ChemicalProps::temperature() const -> real
{
    return T;
}

auto ChemicalProps::pressure() const -> real
{
    return P;
}

auto ChemicalProps::charge() const -> real
{
    const auto Acharge = msystem.formulaMatrixCharge();
    return (Acharge * n.matrix()).sum();
}

auto ChemicalProps::chargeInPhase(StringOrIndex phase) const -> real
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto Az = msystem.formulaMatrixCharge();
    const auto np = n.matrix().segment(offset, length);
    const auto Azp = Az.row(0).segment(offset, length);
    return Azp * np;
}

auto ChemicalProps::chargeAmongSpecies(Indices const& indices) const -> real
{
    const auto Az = msystem.formulaMatrixCharge();
    const auto Azi = Az.row(0)(indices);
    const auto ni = n(indices).matrix();
    return Azi * ni;
}

auto ChemicalProps::elementAmount(StringOrIndex element) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto A = msystem.formulaMatrixElements();
    return A.row(ielement) * n.matrix();
}

auto ChemicalProps::elementAmountInPhase(StringOrIndex element, StringOrIndex phase) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto A = msystem.formulaMatrixElements();
    const auto np = n.matrix().segment(offset, length);
    const auto Aep = A.row(ielement).segment(offset, length);
    return Aep * np;
}

auto ChemicalProps::elementAmountAmongSpecies(StringOrIndex element, Indices const& indices) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto A = msystem.formulaMatrixElements();
    const auto Aei = A.row(ielement)(indices);
    const auto ni = n(indices).matrix();
    return Aei * ni;
}

auto ChemicalProps::elementMass(StringOrIndex element) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmount(ielement);
    return amount * molarmass;
}

auto ChemicalProps::elementMassInPhase(StringOrIndex element, StringOrIndex phase) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmountInPhase(ielement, iphase);
    return amount * molarmass;
}

auto ChemicalProps::elementMassAmongSpecies(StringOrIndex element, Indices const& indices) const -> real
{
    const auto ielement = detail::resolveElementIndexOrRaiseError(msystem, element);
    const auto molarmass = msystem.element(ielement).molarMass();
    const auto amount = elementAmountAmongSpecies(ielement, indices);
    return amount * molarmass;
}

auto ChemicalProps::elementAmounts() const -> ArrayXr
{
    const auto A = msystem.formulaMatrixElements();
    return (A * n.matrix()).array();
}

auto ChemicalProps::elementAmountsInPhase(StringOrIndex phase) const -> ArrayXr
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto A = msystem.formulaMatrixElements();
    const auto Ap = A.middleCols(offset, length);
    const auto np = n.matrix().segment(offset, length);
    return (Ap * np).array();
}

auto ChemicalProps::elementAmountsAmongSpecies(Indices const& indices) const -> ArrayXr
{
    const auto A = msystem.formulaMatrixElements();
    const auto Ai = A(Eigen::all, indices);
    const auto ni = n(indices).matrix();
    return (Ai * ni).array();
}

auto ChemicalProps::componentAmounts() const -> ArrayXr
{
    const auto A = msystem.formulaMatrix();
    return (A * n.matrix()).array();
}

auto ChemicalProps::componentAmountsInPhase(StringOrIndex phase) const -> ArrayXr
{
    const auto iphase = detail::resolvePhaseIndexOrRaiseError(msystem, phase);
    const auto offset = msystem.phases().numSpeciesUntilPhase(iphase);
    const auto length = msystem.phase(iphase).species().size();
    const auto A = msystem.formulaMatrix();
    const auto Ap = A.middleCols(offset, length);
    const auto np = n.matrix().segment(offset, length);
    return (Ap * np).array();
}

auto ChemicalProps::componentAmountsAmongSpecies(Indices const& indices) const -> ArrayXr
{
    const auto A = msystem.formulaMatrix();
    const auto Ai = A(Eigen::all, indices);
    const auto ni = n(indices).matrix();
    return (Ai * ni).array();
}

auto ChemicalProps::speciesAmount(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return n[ispecies];
}

auto ChemicalProps::speciesMass(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return n[ispecies] * msystem.species(ispecies).molarMass();
}

auto ChemicalProps::speciesMoleFraction(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return x[ispecies];
}

auto ChemicalProps::speciesConcentration(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return exp(ln_a[ispecies] - ln_g[ispecies]);
}

auto ChemicalProps::speciesConcentrationLg(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return (ln_a[ispecies] - ln_g[ispecies])/ln10;
}

auto ChemicalProps::speciesConcentrationLn(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return ln_a[ispecies] - ln_g[ispecies];
}

auto ChemicalProps::speciesActivityCoefficient(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return exp(ln_g[ispecies]);
}

auto ChemicalProps::speciesActivityCoefficientLg(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return ln_g[ispecies]/ln10;
}

auto ChemicalProps::speciesActivityCoefficientLn(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return ln_g[ispecies];
}

auto ChemicalProps::speciesActivity(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return exp(ln_a[ispecies]);
}

auto ChemicalProps::speciesActivityLg(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return ln_a[ispecies]/ln10;
}

auto ChemicalProps::speciesActivityLn(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return ln_a[ispecies];
}

auto ChemicalProps::speciesChemicalPotential(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return u[ispecies];
}

auto ChemicalProps::speciesStandardVolume(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return V0[ispecies];
}

auto ChemicalProps::speciesStandardVolumeT(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return VT0[ispecies];
}

auto ChemicalProps::speciesStandardVolumeP(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return VP0[ispecies];
}

auto ChemicalProps::speciesCorrectiveMolarVolume(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return Vxi[ispecies];
}

auto ChemicalProps::speciesStandardGibbsEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return G0[ispecies];
}

auto ChemicalProps::speciesStandardEnthalpy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return H0[ispecies];
}

auto ChemicalProps::speciesStandardEntropy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return (H0[ispecies] - G0[ispecies])/T; // from G0 = H0 - T*S0
}

auto ChemicalProps::speciesStandardInternalEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return H0[ispecies] - P*V0[ispecies]; // from H0 = U0 + P*V0
}

auto ChemicalProps::speciesStandardHelmholtzEnergy(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return G0[ispecies] - P*V0[ispecies]; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ChemicalProps::speciesStandardHeatCapacityConstP(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return Cp0[ispecies];
}

auto ChemicalProps::speciesStandardHeatCapacityConstV(StringOrIndex species) const -> real
{
    const auto ispecies = detail::resolveSpeciesIndexOrRaiseError(msystem, species);
    return VP0[ispecies] == 0.0 ? Cp0[ispecies] : Cp0[ispecies] + T*VT0[ispecies]*VT0[ispecies]/VP0[ispecies]; // from Cv0 = Cp0 + T*VT0*VT0/VP0
}

auto ChemicalProps::speciesAmounts() const -> ArrayXrConstRef
{
    return n;
}

auto ChemicalProps::speciesMasses() const -> ArrayXr
{
    ArrayXr m(n);
    auto i = 0; for(auto const& species : system().species())
        m[i++] *= species.molarMass();
    return m;
}

auto ChemicalProps::speciesMoleFractions() const -> ArrayXrConstRef
{
    return x;
}

auto ChemicalProps::speciesConcentrationsLn() const -> ArrayXr
{
    return ln_a - ln_g;
}

auto ChemicalProps::speciesActivityCoefficientsLn() const -> ArrayXrConstRef
{
    return ln_g;
}

auto ChemicalProps::speciesActivitiesLn() const -> ArrayXrConstRef
{
    return ln_a;
}

auto ChemicalProps::speciesChemicalPotentials() const -> ArrayXrConstRef
{
    return u;
}

auto ChemicalProps::speciesPartialMolarVolumes() const -> ArrayXr
{
    return V0 + Vxi;
}

auto ChemicalProps::speciesStandardVolumes() const -> ArrayXrConstRef
{
    return V0;
}

auto ChemicalProps::speciesStandardVolumesT() const -> ArrayXrConstRef
{
    return VT0;
}

auto ChemicalProps::speciesStandardVolumesP() const -> ArrayXrConstRef
{
    return VP0;
}

auto ChemicalProps::speciesStandardGibbsEnergies() const -> ArrayXrConstRef
{
    return G0;
}

auto ChemicalProps::speciesStandardEnthalpies() const -> ArrayXrConstRef
{
    return H0;
}

auto ChemicalProps::speciesStandardEntropies() const -> ArrayXr
{
    return (H0 - G0)/T; // from G0 = H0 - T*S0
}

auto ChemicalProps::speciesStandardInternalEnergies() const -> ArrayXr
{
    return H0 - P*V0; // from H0 = U0 + P*V0
}

auto ChemicalProps::speciesStandardHelmholtzEnergies() const -> ArrayXr
{
    return G0 - P*V0; // from A0 = U0 - T*S0 = (H0 - P*V0) + (G0 - H0) = G0 - P*V0
}

auto ChemicalProps::speciesStandardHeatCapacitiesConstP() const -> ArrayXrConstRef
{
    return Cp0;
}

auto ChemicalProps::speciesStandardHeatCapacitiesConstV() const -> ArrayXr
{
    return (VP0 == 0.0).select(Cp0, Cp0 + T*VT0*VT0/VP0); // from Cv0 = Cp0 + T*VT0*VT0/VP0
}

auto ChemicalProps::surfaceArea(StringOrIndex const& surface) const -> real
{
    const auto isurface = detail::resolveSurfaceIndexOrRaiseError(msystem, surface);
    const auto numsurfaces = msystem.surfaces().size();
    errorif(isurface >= numsurfaces, "There is no surface in the chemical system with name or index `", stringfy(surface), "`.");
    return msystem.surface(isurface).area(*this);
}

auto ChemicalProps::surfaceAreas() const -> ArrayXr
{
    auto const& surfaces = msystem.surfaces();
    ArrayXr areas(surfaces.size());
    for(auto const& [i, surface] : enumerate(surfaces))
        areas[i] = surface.area(*this);
    return areas;
}

auto ChemicalProps::molarVolume() const -> real
{
    return volume() / amount();
}

auto ChemicalProps::molarVolumeT() const -> real
{
    return volumeT() / amount();
}

auto ChemicalProps::molarVolumeP() const -> real
{
    return volumeP() / amount();
}

auto ChemicalProps::molarGibbsEnergy() const -> real
{
    return gibbsEnergy() / amount();
}

auto ChemicalProps::molarEnthalpy() const -> real
{
    return enthalpy() / amount();
}

auto ChemicalProps::molarEntropy() const -> real
{
    return entropy() / amount();
}

auto ChemicalProps::molarInternalEnergy() const -> real
{
    return internalEnergy() / amount();
}

auto ChemicalProps::molarHelmholtzEnergy() const -> real
{
    return helmholtzEnergy() / amount();
}

auto ChemicalProps::molarHeatCapacityConstP() const -> real
{
    return heatCapacityConstP() / amount();
}

auto ChemicalProps::molarHeatCapacityConstV() const -> real
{
    return heatCapacityConstV() / amount();
}

auto ChemicalProps::specificVolume() const -> real
{
    return volume() / mass();
}

auto ChemicalProps::specificVolumeT() const -> real
{
    return volumeT() / mass();
}

auto ChemicalProps::specificVolumeP() const -> real
{
    return volumeP() / mass();
}

auto ChemicalProps::specificGibbsEnergy() const -> real
{
    return gibbsEnergy() / mass();
}

auto ChemicalProps::specificEnthalpy() const -> real
{
    return enthalpy() / mass();
}

auto ChemicalProps::specificEntropy() const -> real
{
    return entropy() / mass();
}

auto ChemicalProps::specificInternalEnergy() const -> real
{
    return internalEnergy() / mass();
}

auto ChemicalProps::specificHelmholtzEnergy() const -> real
{
    return helmholtzEnergy() / mass();
}

auto ChemicalProps::specificHeatCapacityConstP() const -> real
{
    return heatCapacityConstP() / mass();
}

auto ChemicalProps::specificHeatCapacityConstV() const -> real
{
    return heatCapacityConstV() / mass();
}

auto ChemicalProps::density() const -> real
{
    return mass() / volume();
}

auto ChemicalProps::amount() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).amount(); });
}

auto ChemicalProps::mass() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).mass(); });
}

auto ChemicalProps::volume() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).volume(); });
}

auto ChemicalProps::volumeT() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).volumeT(); });
}

auto ChemicalProps::volumeP() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).volumeP(); });
}

auto ChemicalProps::gibbsEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).gibbsEnergy(); });
}

auto ChemicalProps::enthalpy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).enthalpy(); });
}

auto ChemicalProps::entropy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).entropy(); });
}

auto ChemicalProps::internalEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).internalEnergy(); });
}

auto ChemicalProps::helmholtzEnergy() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).helmholtzEnergy(); });
}

auto ChemicalProps::heatCapacityConstP() const -> real
{
    const auto iend = system().phases().size();
    return Reaktoro::sum(iend, [&](auto i) { return phaseProps(i).heatCapacityConstP(); });
}

auto ChemicalProps::heatCapacityConstV() const -> real
{
    const auto Cp = heatCapacityConstP();
    const auto T = temperature();
    const auto VT = volumeT();
    const auto VP = volumeP();
    return Cp + T*VT*VT/VP;
}

auto ChemicalProps::reactionRate(StringOrIndex reaction) const -> real
{
    const auto ireaction = detail::resolveReactionIndexOrRaiseError(msystem, reaction);
    return msystem.reaction(ireaction).rate(*this);
}

auto ChemicalProps::reactionRates() const -> ArrayXr
{
    auto const& reactions = msystem.reactions();
    ArrayXr rates(reactions.size());
    for(auto const& [i, reaction] : enumerate(reactions))
        rates[i] = reaction.rate(*this);
    return rates;
}

auto ChemicalProps::indicesPhasesWithState(StateOfMatter stateofmatter) const -> Indices
{
    const auto K = msystem.phases().size();
    Indices res;
    res.reserve(K);
    for(auto i = 0; i < K; ++i)
        if(som[i] == stateofmatter)
            res.push_back(i);
    return res;
}

auto ChemicalProps::indicesPhasesWithStates(std::initializer_list<StateOfMatter> stateofmatters) const -> Indices
{
    const auto K = msystem.phases().size();
    Indices res;
    res.reserve(K);
    for(auto i = 0; i < K; ++i)
        if(contains(stateofmatters, som[i]))
            res.push_back(i);
    return res;
}

auto ChemicalProps::indicesPhasesWithFluidState() const -> Indices
{
    return indicesPhasesWithStates({ StateOfMatter::Liquid, StateOfMatter::Gas, StateOfMatter::Supercritical });
}

auto ChemicalProps::indicesPhasesWithSolidState() const -> Indices
{
    return indicesPhasesWithState(StateOfMatter::Solid);
}

auto ChemicalProps::indicesSpeciesInPhasesWithState(StateOfMatter som) const -> Indices
{
    const auto iphases = indicesPhasesWithState(som);
    return msystem.phases().indicesSpeciesInPhases(iphases);
}

auto ChemicalProps::indicesSpeciesInPhasesWithStates(std::initializer_list<StateOfMatter> soms) const -> Indices
{
    const auto iphases = indicesPhasesWithStates(soms);
    return msystem.phases().indicesSpeciesInPhases(iphases);
}

auto ChemicalProps::indicesSpeciesInPhasesWithFluidState() const -> Indices
{
    const auto iphases = indicesPhasesWithFluidState();
    return msystem.phases().indicesSpeciesInPhases(iphases);
}

auto ChemicalProps::indicesSpeciesInPhasesWithSolidState() const -> Indices
{
    const auto iphases = indicesPhasesWithSolidState();
    return msystem.phases().indicesSpeciesInPhases(iphases);
}

auto ChemicalProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto ChemicalProps::output(String const& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

ChemicalProps::operator VectorXr() const
{
    ArrayStream<real> stream;
    serialize(stream);
    return stream.data();
}

ChemicalProps::operator VectorXd() const
{
    ArrayStream<double> stream;
    serialize(stream);
    return stream.data();
}

auto operator<<(std::ostream& out, ChemicalProps const& props) -> std::ostream&
{
    if(props.system().species().empty())
        return out;

    const auto species = props.system().species();
    const auto elements = props.system().elements();
    const auto b   = props.elementAmounts();
    const auto n   = props.speciesAmounts();
    const auto x   = props.speciesMoleFractions();
    const auto lng = props.speciesActivityCoefficientsLn();
    const auto lna = props.speciesActivitiesLn();
    const auto mu  = props.speciesChemicalPotentials();
    const auto G0  = props.speciesStandardGibbsEnergies();
    const auto H0  = props.speciesStandardEnthalpies();
    const auto V0  = props.speciesStandardVolumes();
    const auto S0  = props.speciesStandardEntropies();
    const auto U0  = props.speciesStandardInternalEnergies();
    const auto A0  = props.speciesStandardHelmholtzEnergies();
    const auto Cp0 = props.speciesStandardHeatCapacitiesConstP();
    const auto Cv0 = props.speciesStandardHeatCapacitiesConstV();

    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(props.temperature()), "K" });
    table.add_row({ "Pressure", strfix(props.pressure()*1e-5), "bar" });
    table.add_row({ "Volume", strsci(props.volume()), "m3" });
    table.add_row({ "Gibbs Energy", strfix(props.gibbsEnergy()), "J" });
    table.add_row({ "Enthalpy", strfix(props.enthalpy()), "J" });
    table.add_row({ "Entropy", strfix(props.entropy()), "J/K" });
    table.add_row({ "Internal Energy", strfix(props.internalEnergy()), "J" });
    table.add_row({ "Helmholtz Energy", strfix(props.helmholtzEnergy()), "J" });
    table.add_row({ "Charge", strsci(props.charge()), "mol" });

    table.add_row({ "Element Amount:" }); for(auto i = 0; i < b.size(); ++i) table.add_row({ ":: " + elements[i].symbol(), strsci(b[i]), "mol" });
    table.add_row({ "Species Amount:" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(n[i]), "mol" });
    table.add_row({ "Mole Fraction:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(x[i]), "mol/mol" });
    table.add_row({ "Activity Coefficient:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(exp(lng[i])), "-" });
    table.add_row({ "Activity:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(exp(lna[i])), "-" });
    table.add_row({ "lg(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(lna[i]/ln10), "-" });
    table.add_row({ "ln(Activity):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(lna[i]), "-" });
    table.add_row({ "Chemical Potential:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(mu[i]), "J/mol" });
    table.add_row({ "Standard Volume:", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strsci(V0[i]), "m3/mol" });
    table.add_row({ "Standard Gibbs Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(G0[i]), "J/mol" });
    table.add_row({ "Standard Enthalpy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(H0[i]), "J/mol" });
    table.add_row({ "Standard Entropy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(S0[i]), "J/(mol*K)" });
    table.add_row({ "Standard Internal Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(U0[i]), "J/mol" });
    table.add_row({ "Standard Helmholtz Energy (formation):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(A0[i]), "J/mol" });
    table.add_row({ "Standard Heat Capacity (constant P):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(Cp0[i]), "J/(mol*K)" });
    table.add_row({ "Standard Heat Capacity (constant V):", "", "" }); for(auto i = 0; i < n.size(); ++i) table.add_row({ ":: " + species[i].repr(), strfix(Cv0[i]), "J/(mol*K)" });

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
