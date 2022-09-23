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

#include "Phases.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

GenericPhase::GenericPhase()
{}

GenericPhase::GenericPhase(const StringList& species)
: names(species)
{}

GenericPhase::GenericPhase(const Speciate& elements)
: symbols(elements.symbols)
{}

GenericPhase::GenericPhase(const Speciate& elements, const Exclude& withtags)
: symbols(elements.symbols), excludetags(withtags.tags)
{}

GenericPhase::GenericPhase(const Exclude& withtags)
: excludetags(withtags.tags)
{}

GenericPhase::~GenericPhase()
{}

auto GenericPhase::setName(String name) -> GenericPhase&
{
    phasename = name;
    return *this;
}

auto GenericPhase::setStateOfMatter(StateOfMatter option) -> GenericPhase&
{
    stateofmatter = option;
    return *this;
}

auto GenericPhase::setAggregateState(AggregateState option) -> GenericPhase&
{
    aggregatestate = option;
    return *this;
}

auto GenericPhase::setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GenericPhase&
{
    other_aggregate_states = options;
    return *this;
}

auto GenericPhase::setActivityModel(const ActivityModelGenerator& model) -> GenericPhase&
{
    activity_model = model;
    return *this;
}

auto GenericPhase::setIdealActivityModel(const ActivityModelGenerator& model) -> GenericPhase&
{
    ideal_activity_model = model;
    return *this;
}

auto GenericPhase::named(String name) -> GenericPhase&
{
    return setName(name);
}

auto GenericPhase::set(StateOfMatter option) -> GenericPhase&
{
    return setStateOfMatter(option);
}

auto GenericPhase::set(AggregateState option) -> GenericPhase&
{
    return setAggregateState(option);
}

auto GenericPhase::set(const ActivityModelGenerator& model) -> GenericPhase&
{
    return setActivityModel(model);
}

auto GenericPhase::name() const -> String
{
    return phasename;
}

auto GenericPhase::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GenericPhase::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GenericPhase::additionalAggregateStates() const -> const Vec<AggregateState>&
{
    return other_aggregate_states;
}

auto GenericPhase::species() const -> const Strings&
{
    return names;
}

auto GenericPhase::elements() const -> const Strings&
{
    return symbols;
}

auto GenericPhase::activityModel() const -> const ActivityModelGenerator&
{
    return activity_model;
}

auto GenericPhase::idealActivityModel() const -> const ActivityModelGenerator&
{
    return ideal_activity_model;
}

auto GenericPhase::convert(const Database& db, const Strings& elements) const -> Phase
{
    error(aggregatestate == AggregateState::Undefined,
        "GenericPhase::convert requires an AggregateState value to be specified.\n"
        "Use method GenericPhase::setAggregateState to fix this.");

    auto species = db.speciesWithAggregateState(aggregatestate);

    // If additional aggregate states provided, consider also other species in the database
    for(auto other_aggregate_state : other_aggregate_states)
    {
        auto other_species = db.speciesWithAggregateState(other_aggregate_state);
        if(other_species.size())
            species = concatenate(species, other_species);
    }

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    // Filter out species with provided tags in the exclude function
    if(excludetags.size())
        species = species.withoutTags(excludetags);

    errorif(species.empty(), "Expecting at least one species when defining a phase, but none was provided. Make sure you have listed the species names yourself or used the `speciate` method appropriately.")

    Phase phase;
    phase = phase.withName(phasename);
    phase = phase.withStateOfMatter(stateofmatter);
    phase = phase.withSpecies(species);
    phase = phase.withActivityModel(activity_model(species));
    phase = phase.withIdealActivityModel(ideal_activity_model(species));

    return phase;
}


GenericPhasesGenerator::GenericPhasesGenerator()
{}

GenericPhasesGenerator::GenericPhasesGenerator(const StringList& species)
: names(species)
{}

GenericPhasesGenerator::GenericPhasesGenerator(const Speciate& elements)
: symbols(elements.symbols)
{}

GenericPhasesGenerator::GenericPhasesGenerator(const Speciate& elements, const Exclude& withtags)
: symbols(elements.symbols), excludetags(withtags.tags)
{}

GenericPhasesGenerator::GenericPhasesGenerator(const Exclude& withtags)
: excludetags(withtags.tags)
{}

GenericPhasesGenerator::~GenericPhasesGenerator()
{}

auto GenericPhasesGenerator::setStateOfMatter(StateOfMatter option) -> GenericPhasesGenerator&
{
    stateofmatter = option;
    return *this;
}

auto GenericPhasesGenerator::setAggregateState(AggregateState option) -> GenericPhasesGenerator&
{
    aggregatestate = option;
    return *this;
}

auto GenericPhasesGenerator::setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GenericPhasesGenerator&
{
    other_aggregate_states = options;
    return *this;
}

auto GenericPhasesGenerator::setActivityModel(const ActivityModelGenerator& model) -> GenericPhasesGenerator&
{
    activity_model = model;
    return *this;
}

auto GenericPhasesGenerator::setIdealActivityModel(const ActivityModelGenerator& model) -> GenericPhasesGenerator&
{
    ideal_activity_model = model;
    return *this;
}

auto GenericPhasesGenerator::set(StateOfMatter option) -> GenericPhasesGenerator&
{
    return setStateOfMatter(option);
}

auto GenericPhasesGenerator::set(AggregateState option) -> GenericPhasesGenerator&
{
    return setAggregateState(option);
}

auto GenericPhasesGenerator::set(const ActivityModelGenerator& model) -> GenericPhasesGenerator&
{
    return setActivityModel(model);
}

auto GenericPhasesGenerator::stateOfMatter() const -> StateOfMatter
{
    return stateofmatter;
}

auto GenericPhasesGenerator::aggregateState() const -> AggregateState
{
    return aggregatestate;
}

auto GenericPhasesGenerator::additionalAggregateStates() const -> const Vec<AggregateState>&
{
    return other_aggregate_states;
}

auto GenericPhasesGenerator::species() const -> const Strings&
{
    return names;
}

auto GenericPhasesGenerator::elements() const -> const Strings&
{
    return symbols;
}

auto GenericPhasesGenerator::activityModel() const -> const ActivityModelGenerator&
{
    return activity_model;
}

auto GenericPhasesGenerator::idealActivityModel() const -> const ActivityModelGenerator&
{
    return ideal_activity_model;
}

auto GenericPhasesGenerator::convert(const Database& db, const Strings& elements) const -> Vec<GenericPhase>
{
    error(aggregatestate == AggregateState::Undefined,
        "GenericPhasesGenerator::convert requires an AggregateState value to be specified. "
        "Use method GenericPhasesGenerator::set(AggregateState) to fix this.");

    auto species = db.speciesWithAggregateState(aggregatestate);

    // If additional aggregate states provided, consider also other species in the database
    for(auto other_aggregate_state : other_aggregate_states)
    {
        auto other_species = db.speciesWithAggregateState(other_aggregate_state);
        if(other_species.size())
            species = concatenate(species, other_species);
    }

    species =
        names.size() ? species.withNames(names) :
        symbols.size() ? species.withElements(symbols) :
            species.withElements(elements);

    // Filter out species with provided tags in the exclude function
    if(excludetags.size())
        species = species.withoutTags(excludetags);

    errorif(species.empty(), "Expecting at least one species when defining a list of single-species phases, but none was provided. Make sure you have listed the species names yourself or used the `speciate` method appropriately.")

    Vec<GenericPhase> phases;
    phases.reserve(species.size());
    for(auto&& s : species)
    {
        GenericPhase phase(s.name());
        phase.setName(s.name());
        phase.setStateOfMatter(stateOfMatter());
        phase.setAggregateState(aggregateState());
        phase.setAdditionalAggregateStates(additionalAggregateStates());
        phase.setActivityModel(activityModel());
        phase.setIdealActivityModel(idealActivityModel());

        phases.push_back( phase );
    }

    return phases;
}

Phases::Phases(const Database& db)
: db(db)
{}

auto Phases::add(const GenericPhase& phase) -> void
{
    genericphases.push_back(phase);
}

auto Phases::add(const GenericPhasesGenerator& generator) -> void
{
    generators.push_back(generator);
}

auto Phases::database() const -> const Database&
{
    return db;
}

auto Phases::convert() const -> Vec<Phase>
{
    // Return the element symbols in all stored GenericPhase and GenericPhasesGenerator objects.
    auto collect_all_element_symbols = [&]() -> Strings
    {
        Strings symbols;

        const auto collect_element_symbols_in_genericphase_or_generator = [&](const auto& phase)
        {
            Strings result;
            if(phase.elements().size())
                result = merge(result, phase.elements());
            if(phase.species().size())
                for(auto&& s : db.species().withNames(phase.species()))
                    result = merge(result, s.elements().symbols());
            if(phase.aggregateState() == AggregateState::Aqueous)
                result = merge(result, Strings{"H", "O"}); // ensure both H and O are considered in case there is aqueous phases
            return result;
        };

        for(const auto& phase : genericphases)
            symbols = merge(symbols, collect_element_symbols_in_genericphase_or_generator(phase));

        for(const auto& generator : generators)
            symbols = merge(symbols, collect_element_symbols_in_genericphase_or_generator(generator));

        return symbols;
    };

    // Return all given GenericPhase objects together with the generated ones.
    auto collect_all_generic_phases = [&](const Strings& symbols) -> Vec<GenericPhase>
    {
        Vec<GenericPhase> collected(genericphases);

        for(auto generator : generators)
        {
            const Vec<GenericPhase> generated = generator.convert(db, symbols);
            collected.insert(collected.end(), generated.begin(), generated.end());
        }

        return collected;
    };

    // Replace duplicate phase names with unique names.
    auto fix_duplicate_phase_names = [](Vec<GenericPhase>& phases)
    {
        Strings phasenames = vectorize(phases, RKT_LAMBDA(x, x.name()));
        phasenames = makeunique(phasenames, "!");
        auto i = 0;
        for(auto& phase : phases)
            phase.setName(phasenames[i++]);
    };

    Strings symbols = collect_all_element_symbols();

    Vec<GenericPhase> allgenericphases = collect_all_generic_phases(symbols);

    fix_duplicate_phase_names(allgenericphases);

    Vec<Phase> phases;
    phases.reserve(allgenericphases.size());
    for(const auto& genericphase : allgenericphases)
        phases.push_back(genericphase.convert(db, symbols));

    return phases;
}

Phases::operator Vec<Phase>() const
{
    return convert();
}

} // namespace Reaktoro
