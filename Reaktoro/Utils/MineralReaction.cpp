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

#include "MineralReaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

namespace Reaktoro {
namespace detail {

// /// Convert a vector of MineralReactionRateModel objects to a vector of ReactionRateModel objects.
// auto convert(Strings const& minerals, Vec<MineralReactionRateModel> const& models) -> Vec<ReactionRateModel>
// {
//     errorif(minerals.size() != models.size(), "Expecting same number of minerals and mineral reaction rate models.");

//     Vec<ReactionRateModel> ratemodels(models.size());

//     for(auto i = 0; i < models.size(); ++i)
//     {
//         ratemodels[i] = [=](ChemicalProps const& props) -> ReactionRate
//         {
//             auto const& system = props.system();
//             auto const& aprops = getAqueousProps(system); // get the AqueousProps object corresponding to the ChemicalSystem object in `props`, which will be reused by all mineral reactions!

//             thread_local auto const iminerals = vectorize(minerals, RKT_LAMBDA(x, aprops.saturationSpecies().indexWithName(x)));
//             thread_local auto const imineralphases = vectorize(minerals, RKT_LAMBDA(x, system.phases().indexWithName(x)));
//             thread_local auto const imineralsurfaces = vectorize(minerals, RKT_LAMBDA(x, system.surfaces().indexWithName(x)));

//             auto const& T = props.temperature();
//             auto const& P = props.pressure();
//             auto const& pH = aprops.pH();
//             auto const& Omega = aprops.saturationRatio(iminerals[i]);
//             auto const& area = props.surfaceArea(imineralsurfaces[i]);

//             const auto args = MineralReactionRateModelArgs{ props, aprops, T, P, pH, Omega, area };

//             // Evaluate the mineral reaction rate model and switch sign because
//             // Reaktoro's convention for reaction rate is positive when the
//             // reaction proceeds from left to right (and this is how the mineral
//             // reaction is represented, with the mineral on the left side, its
//             // dissolution rate, from left to right, should be positive).
//             return -models[i](args);
//         };
//     }

//     // Ensure the first rate model updates AqueousProps object before it is used by every mineral rate model!
//     ratemodels[0] = [=](ChemicalProps const& props) mutable -> ReactionRate
//     {
//         getAqueousProps(props.system()).update(props); // update the AqueousProps object corresponding to the ChemicalSystem object in `props`; note this is done once for the first mineral rate model, and updated state reused by all other minerals for performance reasons!
//         return ratemodels[0](props);
//     };

//     return ratemodels;
// }

/// Convert a MineralReactionRateModel object to a ReactionRateModel object.
auto convert(String const& mineral, MineralReactionRateModel const& model) -> ReactionRateModel
{
    errorif(!model, "Expecting an initialized MineralReactionRateModel object when converting it to a ReactionRateModel");

    return [=](ChemicalProps const& props) -> ReactionRate
    {
        auto const& aprops = AqueousProps::compute(props);
        auto const& T = props.temperature();
        auto const& P = props.pressure();
        auto const& pH = aprops.pH();
        auto const& Omega = aprops.saturationRatio(mineral);
        auto const& area = props.surfaceArea(mineral);

        const auto args = MineralReactionRateModelArgs{ props, aprops, T, P, pH, Omega, area };

        // Evaluate the mineral reaction rate model and switch sign because
        // Reaktoro's convention for reaction rate is positive when the
        // reaction proceeds from left to right (and this is how the mineral
        // reaction is represented, with the mineral on the left side, its
        // dissolution rate, from left to right, should be positive).
        return -model(args);
    };
}

} // namespace detail

MineralReaction::MineralReaction(String const& mineral)
: GeneralReaction(mineral)
{}

auto MineralReaction::setRateModel(MineralReactionRateModel const& model) -> MineralReaction&
{
    ReactionRateModel converted = detail::convert(mineral(), model);
    GeneralReaction::setRateModel(converted);
    return *this;
}

auto MineralReaction::setRateModel(MineralReactionRateModelGenerator const& model_generator) -> MineralReaction&
{
    const auto mineralname = mineral();
    ReactionRateModelGenerator converted = [=](ReactionRateModelGeneratorArgs args)
    {
        MineralReactionRateModel model = model_generator(args);
        return detail::convert(mineralname, model);
    };
    GeneralReaction::setRateModel(converted);
    return *this;
}

auto MineralReaction::mineral() const -> String const&
{
    return GeneralReaction::name();
}
// MineralReaction::MineralReaction(String const& mineral)
// : m_minerals(minerals), m_mineral_rate_model_generators(m_minerals.size())
// {
// }

// auto MineralReaction::setRateModel(MineralReactionRateModelGenerator const& generator) -> MineralReaction&
// {
//     const auto size = m_mineral_rate_model_generators.size();
//     m_mineral_rate_model_generators.assign(size, generator);
//     return *this;
// }

// auto MineralReaction::setRateModel(String const& mineral, MineralReactionRateModelGenerator const& generator) -> MineralReaction&
// {
//     errorif(!generator, "You are trying to specify a non-initialized mineral reaction rate generator.");
//     const auto idx = index(m_minerals, mineral);
//     errorif(idx >= m_minerals.size(), "You did not specify mineral `", mineral, "` in the list of minerals when creating the MineralReaction object, e.g., `MineralReaction(\"Calcite Dolomite Quartz\")`.");
//     m_mineral_rate_model_generators[idx] = generator;
//     return *this;
// }

// auto MineralReaction::operator()(SpeciesList const& species) const -> Vec<Reaction>
// {
//     for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
//     {
//         const auto imineral = species.findWithName(m_minerals[i]);
//         errorif(!generator, "You forgot to set a mineral reaction rate model for mineral `", m_minerals[i], "` and maybe for all other minerals as well. Use method MineralReaction::setRateModel to fix this.");
//         errorif(imineral >= species.size(), "There is no mineral with name `", m_minerals[i], "` in the chemical system.");
//     }

//     const auto num_minerals = m_minerals.size();

//     Vec<MineralReactionRateModel> mineral_reaction_rate_models(num_minerals);
//     for(auto&& [i, generator] : enumerate(m_mineral_rate_model_generators))
//         mineral_reaction_rate_models[i] = generator(m_minerals[i], species);

//     const auto reaction_rate_models = detail::convert(m_minerals, mineral_reaction_rate_models);

//     Vec<Reaction> reactions(num_minerals);
//     for(auto i = 0; i < num_minerals; ++i)
//     {
//         const auto imineral = species.index(m_minerals[i]);
//         const auto mineralspecies = species[imineral];
//         reactions[i] = Reaction()
//             .withName(m_minerals[i])
//             .withEquation({{ {mineralspecies, -1.0} }})
//             .withRateModel(reaction_rate_models[i])
//             ;
//     }

//     return reactions;
// }

} // namespace Reaktoro