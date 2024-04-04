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

#include "AqueousProps.hpp"

// C++ includes
#include <fstream>

// cpp-tabulate includes
#include <tabulate/table.hpp>
using namespace tabulate;

// Optima includes
#include <Optima/Echelonizer.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Warnings.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Utils.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterUtils.hpp>

namespace Reaktoro {
namespace {

/// Return the index of the first aqueous phase in the system.
auto indexAqueousPhase(ChemicalSystem const& system) -> Index
{
    const auto aqueous_phases = system.phases().withAggregateState(AggregateState::Aqueous);
    warning(aqueous_phases.size() > 1,
        "While creating an AqueousProps object, it has been detected ",
        "more than one aqueous phase in the system. The AqueousProps object "
        "created will correspond to the first aqueous phase found.");
    const auto idx = system.phases().findWithAggregateState(AggregateState::Aqueous);
    error(idx >= system.phases().size(),
        "Could not create an AqueousProps object because there is no "
        "phase in the system with aggregate state value AggregateState::Aqueous.");
    return idx;
}

// Return a chemical potential function for a non-aqueous species using a given
// activity model. Note: this method is relevant for computation of saturation indices.
auto chemicalPotentialModel(Species const& species, ActivityModelGenerator const& generator) -> Fn<real(ChemicalProps const&)>
{
    const auto activitymodel = generator({species}); // TODO: Use .withMemoization() here to avoid full recomputations when same T and P are given.
    const auto R = universalGasConstant;
    const auto x = ArrayXr{{1.0}}; // the mole fraction of the single species in a pure phase
    auto actprops = ActivityProps::create(1);

    return [=](ChemicalProps const& props) mutable -> real
    {
        const auto T = props.temperature();
        const auto P = props.pressure();
        activitymodel(actprops, {T, P, x}); // evaluate the activity model
        const auto G0 = species.standardThermoProps(T, P).G0;
        const auto ln_a = actprops.ln_a[0];
        return G0 + R*T*ln_a;
    };
}

// Return a default chemical potential function for a non-aqueous species. If
// the species exists in the given chemical system, then the chemical potential
// model reuses the available chemical potential in the `ChemicalProps` object
// provided as an argument to the model function. If not, the standard Gibbs
// energy `G0` of the species is computed and used as chemical potential. If
// the species is a gas, then `G0 + RT*ln(Pbar)` is used instead. This default
// behavior implies that the species constitute a pure ideal gas or solid
// phase. Note: this method is relevant for computation of saturation indices.
auto defaultChemicalPotentialModel(Species const& species, ChemicalSystem const& system) -> Fn<real(ChemicalProps const&)>
{
    const auto numspecies = system.species().size();
    const auto ispecies = system.species().find(species.name());

    // Case I: when species exists in the chemical system
    if(ispecies < numspecies)
        return [=](ChemicalProps const& props) -> real
        {
            return props.speciesChemicalPotential(ispecies);
        };
    // Case II: when species is a gas and it does not exist in the chemical system
    else if(species.aggregateState() == AggregateState::Gas)
        return [=](ChemicalProps const& props) -> real
        {
            const auto T = props.temperature();
            const auto P = props.pressure();
            const auto Pbar = P*1e-5; // from Pa to bar
            const auto RT = universalGasConstant * T;
            const auto G0 = species.standardThermoProps(T, P).G0;
            return G0 + RT*log(Pbar);
        };
    // Case III: when species is not a gas and it does not exist in the chemical system
    else
        return [=](ChemicalProps const& props) -> real
        {
            const auto T = props.temperature();
            const auto P = props.pressure();
            const auto G0 = species.standardThermoProps(T, P).G0;
            return G0;
        };
}

// Return a vector with default chemical potential functions for every given chemical species.
auto defaultChemicalPotentialModels(SpeciesList const& nonaqueous, ChemicalSystem const& system) -> Vec<Fn<real(ChemicalProps const&)>>
{
    return vectorize(nonaqueous, RKT_LAMBDA(x, defaultChemicalPotentialModel(x, system)));
}

// Return the index of species H+ if found, otherwise the index of H3O+.
auto findHydrogenIon(SpeciesList const& species) -> Index
{
    const auto idx = species.findWithFormula("H+");
    if(idx < species.size())
        return idx;
    else return species.findWithFormula("H3O+");
}

} // namespace

struct AqueousProps::Impl
{
    /// The chemical system in which the aqueous phase is.
    const ChemicalSystem system;

    /// The index of the underlying Phase object for the aqueous phase in the system.
    const Index iphase;

    /// The underlying Phase object for the aqueous phase in the system.
    const Phase phase;

    /// The phase as an aqueous solution.
    const AqueousMixture aqsolution;

    /// The index of the aqueous solvent species H2O in the aqueous phase (not in the system!)
    const Index iH2O;

    /// The index of chemical species H+ or H3O+ in the aqueous phase (not in the system!)
    const Index iH;

    /// The chemical properties of the system.
    ChemicalProps props;

    /// The state of the aqueous solution.
    AqueousMixtureState aqstate;

    /// The amounts of the species in the aqueous phase (to be used with echelonizer - not for any computation, since it does not have autodiff propagation!).
    VectorXd naq;

    /// The chemical potentials of the elements in the aqueous phase
    VectorXr lambda;

    /// The non-aqueous species in the database for which saturation indices are calculated.
    SpeciesList nonaqueous;

    /// The formula matrix of the aqueous species in the aqueous phase.
    MatrixXd Aaqs;

    /// The formula matrix of the non-aqueous species for the computation of saturation indices.
    MatrixXd Anon;

    /// The echelon form of the formula matrix `Aaqs` of the aqueous species.
    Optima::Echelonizer echelonizer;

    /// The chemical potential models for the non-aqueous species (as if they were pure phases) for the computation of their saturation indices.
    Vec<Fn<real(ChemicalProps const&)>> chemical_potential_models;

    /// The alkalinity contribution factors of some aqueous species based on the alkalinity model of Wolf-Gladrow et al. (2007)
    Pairs<double, Index> alkalinity_factors;

    Impl(ChemicalSystem const& system)
    : system(system),
      iphase(indexAqueousPhase(system)),
      phase(system.phase(iphase)),
      props(system),
      aqsolution(phase.species()),
      iH2O(phase.species().findWithFormula("H2O")),
      iH(findHydrogenIon(phase.species()))
    {
        const auto Naq = phase.species().size();

        error(iH2O >= Naq, "Cannot create AqueousProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H2O.");

        error(iH >= Naq, "Cannot create AqueousProps object for phase ", phase.name(), " "
            "because it does not contain a species with formula H+ or H3O+.");

        // The symbols of the elements in the aqueous phase
        const auto symbols = vectorize(phase.elements(), RKT_LAMBDA(x, x.symbol()));

        // The aqueous species in the aqueous phase
        auto const& aqspecies = phase.species();

        // Collect the species from the database that contains the elements in the aqueous phase
        const auto species_same_elements = system.database().species().withElements(symbols);

        // Collect the non-aqueous species from the database that contains the elements in the aqueous phase
        nonaqueous = removefn(species_same_elements, RKT_LAMBDA(x, x.aggregateState() == AggregateState::Aqueous));

        // Ensure non-aqueous species are sorted by aggregate state (gases, solids, etc)
        std::sort(nonaqueous.begin(), nonaqueous.end(),
            [](auto l, auto r)
                { return l.aggregateState() < r.aggregateState(); });

        // Assemble the formula matrices of the aqueous and non-aqueous species w.r.t. elements in the aqueous phase
        Aaqs = detail::assembleFormulaMatrix(phase.species(), phase.elements());
        Anon = detail::assembleFormulaMatrix(nonaqueous, phase.elements());

        // Initialize the chemical potential models for the non-aqueous species, as if they were pure phases
        chemical_potential_models = defaultChemicalPotentialModels(nonaqueous, system);

        // Initialize the aqueous state properties
        aqstate.T = NaN;
        aqstate.P = NaN;
        aqstate.rho = NaN;
        aqstate.epsilon = NaN;
        aqstate.Ie = NaN;
        aqstate.Is = NaN;
        aqstate.m.setConstant(Naq, NaN);
        aqstate.ms.setConstant(Naq, NaN);

        // Compute the initial echelon form of formula matrix `Aaqs`
        echelonizer.compute(Aaqs);

        // Initialize the contribution factors of some aqueous species towards the alkalinity of the solution (Wolf-Gladrow et al., 2007)
        alkalinity_factors = {
            { 1.0, aqspecies.findWithFormula("Na+")   },
            { 2.0, aqspecies.findWithFormula("Mg+2")  },
            { 2.0, aqspecies.findWithFormula("Ca+2")  },
            { 1.0, aqspecies.findWithFormula("K+")    },
            { 2.0, aqspecies.findWithFormula("Sr+2")  },
            {-1.0, aqspecies.findWithFormula("Cl-")   },
            {-1.0, aqspecies.findWithFormula("Br-")   },
            {-1.0, aqspecies.findWithFormula("NO3-")  },
            {-1.0, aqspecies.findWithFormula("H3PO4") }, // TPO4 in reference.
            {-1.0, aqspecies.findWithFormula("H2PO4-")}, // TPO4 in reference.
            {-1.0, aqspecies.findWithFormula("HPO4-2")}, // TPO4 in reference. Note that factor here is -1 as in the paper.
            {-1.0, aqspecies.findWithFormula("PO4-3") }, // TPO4 in reference. Note that factor here is -1 as in the paper.
            { 1.0, aqspecies.findWithFormula("NH3")   }, // TNH3 in reference.
            { 1.0, aqspecies.findWithFormula("NH4+")  }, // TNH3 in reference.
            {-2.0, aqspecies.findWithFormula("SO4-2") }, // TSO4 in reference.
            {-2.0, aqspecies.findWithFormula("HSO4-") }, // TSO4 in reference. Note that factor here is -2 as in the paper.
            {-1.0, aqspecies.findWithFormula("F-")    }, // THF in reference.
            {-1.0, aqspecies.findWithFormula("HF")    }, // THF in reference.
            {-1.0, aqspecies.findWithFormula("NO2-")  }, // THNO2 in reference.
            {-1.0, aqspecies.findWithFormula("HNO2")  }, // THNO2 in reference. Note that factor here is -1 as in the paper.
        };
    }

    Impl(ChemicalState const& state)
    : Impl(state.system())
    {
        update(state);
    }

    Impl(ChemicalProps const& props)
    : Impl(props.system())
    {
        update(props);
    }

    auto setActivityModel(StringOrIndex const& species, ActivityModelGenerator const& generator) -> void
    {
        const auto i = detail::resolveSpeciesIndex(nonaqueous, species);
        errorif(i >= nonaqueous.size(), "Could not set the activity model of species with name or index `", stringfy(species), "`. "
            "This species must be non-aqueous and exist in the thermodynamic database. It must also be composed of chemical elements "
            "present in the aqueous phase. This error will occur, for example, if you are calculating the saturation ratio of Quartz (SiO2) "
            "but the aqueous phase has no species with element Si.");
        chemical_potential_models[i] = chemicalPotentialModel(nonaqueous[i], generator);
    }

    auto update(ChemicalState const& state) -> void
    {
        props.update(state);
        update(props);
    }

    auto update(ChemicalProps const& cprops) -> void
    {
        // Auxiliary variables
        auto const& aqprops = cprops.phaseProps(iphase);
        auto const& T = aqprops.temperature();
        auto const& P = aqprops.pressure();
        auto const& x = aqprops.speciesMoleFractions();

        // Update the internal properties of the chemical system
        props = cprops;

        // Update the internal aqueous state object
        aqstate = aqsolution.state(T, P, x);

        // Update auxiliary vector naq to be used in the echelonization below
        naq = aqprops.speciesAmounts();

        // Update the echelon form and also the list of basic species
        echelonizer.updateWithPriorityWeights(naq);

        // Compute chemical potentials of the elements in the aqueous phase
        const auto u = aqprops.speciesChemicalPotentials();
        const auto ib = echelonizer.indicesBasicVariables();
        const auto R = echelonizer.R();
        const auto Rb = R.topRows(ib.size());
        const VectorXr ub = u(ib);
        lambda = Rb.transpose() * ub;
    }

    auto temperature() const -> real
    {
        return props.temperature();
    }

    auto pressure() const -> real
    {
        return props.pressure();
    }

    auto waterAmount() const -> real
    {
        return props.speciesAmount(iH2O);
    }

    auto waterMass() const -> real
    {
        return props.speciesMass(iH2O);
    }

    auto charge() const -> real
    {
        return props.chargeInPhase(iphase);
    }

    auto chargeMolality() const -> real
    {
        return charge() / waterMass();
    }

    auto elementMolality(StringOrIndex const& symbol) const -> real
    {
        const auto idx = detail::resolveElementIndexOrRaiseError(phase, symbol);
        auto const& m = aqstate.m.matrix();
        return Aaqs.row(idx) * m;
    }

    auto elementMolalities() const -> ArrayXr
    {
        const auto E = phase.elements().size();
        auto const& m = aqstate.m.matrix();
        return Aaqs.topRows(E) * m;
    }

    auto speciesMolality(StringOrIndex const& name) const -> real
    {
        const auto idx = detail::resolveSpeciesIndexOrRaiseError(phase, name);
        return aqstate.m[idx];
    }

    auto speciesMolalities() const -> ArrayXr
    {
        return aqstate.m;
    }

    auto ionicStrength() const -> real
    {
        return aqstate.Ie;
    }

    auto ionicStrengthStoichiometric() const -> real
    {
        return aqstate.Is;
    }

    auto pH() const -> real
    {
        auto const& aqprops = props.phaseProps(iphase);
        auto const& ln_aH = aqprops.speciesActivitiesLn()[iH];
        return -ln_aH/ln10;
    }

    auto pE() const -> real
    {
        const auto T = props.temperature();
        const auto E = phase.elements().size();
        const auto lambdaZ = lambda[E];
        const auto RT = universalGasConstant * T;
        const auto res = lambdaZ/(RT*ln10);
        return res;
    }

    auto Eh() const -> real
    {
        const auto T = props.temperature();
        const auto RT = universalGasConstant * T;
        const auto F = faradayConstant;
        const auto res = ln10*RT/F*pE();
        return res;
    }

    auto alkalinity() const -> real
    {
        auto const& aqspecies = phase.species();
        auto const& aqprops = props.phaseProps(iphase);
        auto const& naq = aqprops.speciesAmounts();
        auto const num_aqspecies = aqspecies.size();

        real alkalinity = 0.0;
        for(const auto& [factor, iaqspecies] : alkalinity_factors)
            if(iaqspecies < num_aqspecies)
                alkalinity += factor * naq[iaqspecies];

        // Convert to eq/L
        auto const m3_to_liter = 1000.0;
        auto V = aqprops.volume();  // in m3
        auto V_in_liter = V * m3_to_liter;

        if(V <= 0.0)
        {
            warningif(V <= 0.0 && Warnings::isEnabled(548),
                "The chosen thermodynamic model for the aqueous phase produced volume ", V_in_liter, " L, which is not physical. "
                "As an approximation, alkalinity will be computed using density of pure water provided by Wagner & Pruss (2002) model as the density of the solution. "
                "Was the Extended UNIQUAC the chosen model? It currently does not compute densities for lack of standard molar volumes of the species. "
                "Disable this warning message with Warnings.disable(345) in Python and Warnings::disable(345) in C++.");

            auto const T = temperature();
            auto const P = pressure();
            auto const mass = props.phaseProps(iphase).mass(); // aqueous solution mass in kg
            auto const rho = waterDensityWagnerPruss(T, P, StateOfMatter::Liquid); // density in kg/m3
            V = mass / rho;
            V_in_liter = V * m3_to_liter;
        }

        alkalinity /= V_in_liter;

        return alkalinity;
    }

    auto saturationRatioLn(StringOrIndex const& species) const -> real
    {
        const auto i = detail::resolveSpeciesIndex(nonaqueous, species);
        errorif(i >= nonaqueous.size(), "It was not possible to calculate the saturation ratio of "
            "species with name or index `", stringfy(species), "`. This species must be non-aqueous "
            "and exist in the thermodynamic database. It must also be composed of chemical elements "
            "present in the aqueous phase. This error will occur, for example, if you are calculating "
            "the saturation ratio of Quartz (SiO2) but the aqueous phase has no species with element Si.");
        const auto RT = universalGasConstant * props.temperature();
        const auto ui = chemical_potential_models[i](props);
        const auto li = Anon.col(i).dot(lambda);
        const auto lnOmegai = (li - ui)/RT;
        return lnOmegai;
    }

    auto saturationRatiosLn() const -> ArrayXr
    {
        const auto RT = universalGasConstant * props.temperature();
        const auto num_nonaqueous = nonaqueous.size();
        ArrayXr lnOmega(num_nonaqueous);
        lnOmega = Anon.transpose() * lambda;
        for(auto i = 0; i < num_nonaqueous; ++i)
            lnOmega[i] -= chemical_potential_models[i](props);
        lnOmega /= RT;
        return lnOmega;
    }
};

AqueousProps::AqueousProps(ChemicalSystem const& system)
: pimpl(new Impl(system))
{}

AqueousProps::AqueousProps(ChemicalState const& state)
: pimpl(new Impl(state))
{}

AqueousProps::AqueousProps(ChemicalProps const& props)
: pimpl(new Impl(props))
{}

AqueousProps::AqueousProps(AqueousProps const& other)
: pimpl(new Impl(*other.pimpl))
{}

AqueousProps::~AqueousProps()
{}

auto AqueousProps::operator=(AqueousProps other) -> AqueousProps&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto AqueousProps::compute(ChemicalProps const& props) -> AqueousProps const&
{
    using SystemID = Index;
    using CacheEntry = Pair<AqueousProps, const ChemicalProps*>;
    thread_local Map<SystemID, CacheEntry> cache;

    const auto systemid = props.system().id();
    if(auto it = cache.find(systemid); it != cache.end())
    {
        auto& [aprops, pprops] = it->second;

        auto new_props_stateid = props.stateid() != aprops.props().stateid();
        auto new_props_object = &props != pprops;

        if(new_props_stateid || new_props_object) // not enough to check only stateid; two different ChemicalProps with same stateid but different state would result in the cached aprops to be returned.
        {
            aprops.update(props);
            pprops = &props;
        }

        return aprops;
    }

    const auto [it, _] = cache.emplace(systemid, CacheEntry{ AqueousProps(props), &props });
    return it->second.first;
}

auto AqueousProps::setActivityModel(StringOrIndex const& species, ActivityModelGenerator const& generator) -> void
{
    pimpl->setActivityModel(species, generator);
}

auto AqueousProps::update(ChemicalState const& state) -> void
{
    pimpl->update(state);
}

auto AqueousProps::update(ChemicalProps const& props) -> void
{
    pimpl->update(props);
}

auto AqueousProps::temperature() const -> real
{
    return pimpl->temperature();
}

auto AqueousProps::pressure() const -> real
{
    return pimpl->pressure();
}

auto AqueousProps::waterAmount() const -> real
{
    return pimpl->waterAmount();
}

auto AqueousProps::waterMass() const -> real
{
    return pimpl->waterMass();
}

auto AqueousProps::charge() const -> real
{
    return pimpl->charge();
}

auto AqueousProps::chargeMolality() const -> real
{
    return pimpl->chargeMolality();
}

auto AqueousProps::elementMolality(StringOrIndex const& symbol) const -> real
{
    return pimpl->elementMolality(symbol);
}

auto AqueousProps::elementMolalities() const -> ArrayXr
{
    return pimpl->elementMolalities();
}

auto AqueousProps::speciesMolality(StringOrIndex const& name) const -> real
{
    return pimpl->speciesMolality(name);
}

auto AqueousProps::speciesMolalities() const -> ArrayXr
{
    return pimpl->speciesMolalities();
}

auto AqueousProps::ionicStrength() const -> real
{
    return pimpl->ionicStrength();
}

auto AqueousProps::ionicStrengthEffective() const -> real
{
    return pimpl->ionicStrength();
}

auto AqueousProps::ionicStrengthStoichiometric() const -> real
{
    return pimpl->ionicStrengthStoichiometric();
}

auto AqueousProps::pH() const -> real
{
    return pimpl->pH();
}

auto AqueousProps::pE() const -> real
{
    return pimpl->pE();
}

auto AqueousProps::Eh() const -> real
{
    return pimpl->Eh();
}

auto AqueousProps::alkalinity() const -> real
{
    return pimpl->alkalinity();
}

auto AqueousProps::saturationSpecies() const -> SpeciesList
{
    return pimpl->nonaqueous;
}

auto AqueousProps::saturationIndex(StringOrIndex const& species) const -> real
{
    return pimpl->saturationRatioLn(species) / ln10;
}

auto AqueousProps::saturationIndices() const -> ArrayXr
{
    return pimpl->saturationRatiosLn() / ln10; // in log10 scale
}

auto AqueousProps::saturationRatio(StringOrIndex const& species) const -> real
{
    return exp(pimpl->saturationRatioLn(species));
}

auto AqueousProps::saturationRatios() const -> ArrayXr
{
    return pimpl->saturationRatiosLn().exp();
}

auto AqueousProps::saturationRatiosLn() const -> ArrayXr
{
    return pimpl->saturationRatiosLn();
}

auto AqueousProps::props() const -> ChemicalProps const&
{
    return pimpl->props;
}

auto AqueousProps::system() const -> ChemicalSystem const &
{
    return pimpl->system;
}

auto AqueousProps::phase() const -> Phase const&
{
    return pimpl->phase;
}

auto AqueousProps::output(std::ostream& out) const -> void
{
    out << *this;
}

auto AqueousProps::output(String const& filename) const -> void
{
    auto out = std::ofstream(filename);
    out << *this;
}

auto AqueousProps::saturationIndexLn(StringOrIndex const& species) const -> real
{
    errorif(true, "Method AqueousProps::saturationIndexLn has been deprecated. Rely on the use of saturationIndex(species) instead.");
    return {};
}

auto AqueousProps::saturationIndexLg(StringOrIndex const& species) const -> real
{
    errorif(true, "Method AqueousProps::saturationIndexLg has been deprecated. Rely on the use of saturationIndex(species) instead.");
    return {};
}

auto AqueousProps::saturationIndicesLn() const -> ArrayXr
{
    errorif(true, "Method AqueousProps::saturationIndicesLn has been deprecated. Rely on the use of saturationIndices() instead.");
    return {};
}

auto AqueousProps::saturationIndicesLg() const -> ArrayXr
{
    errorif(true, "Method AqueousProps::saturationIndicesLg has been deprecated. Rely on the use of saturationIndices() instead.");
    return {};
}

auto operator<<(std::ostream& out, AqueousProps const& props) -> std::ostream&
{
    const auto elements = props.phase().elements();
    const auto species = props.phase().species();
    const auto ms = props.speciesMolalities();
    const auto me = props.elementMolalities();
    const auto lgOmega = props.saturationIndices();
    assert(species.size() == ms.size());
    assert(elements.size() == me.size());
    Table table;
    table.add_row({ "Property", "Value", "Unit" });
    table.add_row({ "Temperature", strfix(props.temperature()), "K" });
    table.add_row({ "Pressure", strfix(props.pressure()*1e-5), "bar" });
    table.add_row({ "Ionic Strength (Effective)", strfix(props.ionicStrength()), "molal" });
    table.add_row({ "Ionic Strength (Stoichiometric)", strfix(props.ionicStrengthStoichiometric()), "molal" });
    table.add_row({ "pH", strfix(props.pH()), "" });
    table.add_row({ "pE", strfix(props.pE()), "" });
    table.add_row({ "Eh", strfix(props.Eh()), "V" });
    table.add_row({ "Alkalinity", strfix(props.alkalinity()), "eq/L" });
    table.add_row({ "Charge Molality", strsci(props.chargeMolality()), "molal" });
    table.add_row({ "Element Molality:" });
    for(auto i = 0; i < elements.size(); ++i)
        if(elements[i].symbol() != "H" && elements[i].symbol() != "O")
            table.add_row({ ":: " + elements[i].symbol(), strsci(me[i]), "molal" });
    table.add_row({ "Species Molality:" });
    for(auto i = 0; i < species.size(); ++i)
        if(species[i].formula().str() != "H2O")
            table.add_row({ ":: " + species[i].repr(), strsci(ms[i]), "molal" });
    table.add_row({ "Saturation Indices:" });
    for(auto [i, species] : enumerate(props.saturationSpecies()))
        table.add_row({ ":: " + species.repr(), strfix(((lgOmega[i] + 1000)) - 1000), "-" }); // + 1000 - 1000 as a trick to transform -1e15 into 0.0
    table.add_row({ "Saturation Ratios:" });
    for(auto [i, species] : enumerate(props.saturationSpecies()))
        table.add_row({ ":: " + species.repr(), strsci(pow(10.0, lgOmega[i])), "-" });

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
