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

#include "DatabaseParser.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Core/Data.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/StandardThermoModel.hpp>
#include <Reaktoro/Serialization.hpp>

namespace Reaktoro {

struct DatabaseParser::Impl
{
    ///< The Species objects in the database.
    SpeciesList species_list;

    ///< The Element objects in the database.
    ElementList element_list;

    ///< The database contents parsed from YAML or JSON into a Data object.
    Data doc;

    /// Construct a default DatabaseParser::Impl object.
    Impl()
    {}

    /// Construct a DatabaseParser::Impl object with given Data object.
    Impl(const Data& doc)
    : Impl(doc, {})
    {}

    /// Construct a DatabaseParser::Impl object with given Data object.
    Impl(const Data& doc, const ElementList& elements)
    : doc(doc), element_list(elements)
    {
        errorif(!doc.isDict(), "Could not understand your YAML or JSON database file with content:\n", doc.repr(), "\n",
            "Repeating the error message here in case the above printed content is too long.\n",
            "Could not understand your YAML or JSON database file with the above content.\n",
            "Are you forgetting to add the list of chemical species inside a Species YAML or JSON map?\n",
            "Please check other Reaktoro's YAML or JSON databases to identify what is not conforming.");

        // if(doc.exists("Elements"))
        //     for(auto const& child : doc["Elements"].asDict())
        //         addElement(child.first, child.second);

        // if(doc.exists("Species"))
        //     for(auto const& child : doc["Species"].asDict())
        //         addSpecies(child.first, child.second);

        if(doc.exists("Elements"))
        {
            if(doc["Elements"].isDict())
                for(auto const& child : doc["Elements"].asDict())
                    addElement(child.first, child.second);
            else if(doc["Elements"].isList())
                for(auto const& child : doc["Elements"].asList())
                    addElement(child["Symbol"].asString(), child);
            else errorif(true, "Expecting the `Elements` section in your YAML or JSON database to be either a list or dictionary. Please check other Reaktoro databases in either YAML or JSON format and replicate the structure.");
        }

        if(doc.exists("Species"))
        {
            if(doc["Species"].isDict())
                for(auto const& child : doc["Species"].asDict())
                    addSpecies(child.first, child.second);
            else if(doc["Species"].isList())
                for(auto const& child : doc["Species"].asList())
                    addSpecies(child["Name"].asString(), child);
            else errorif(true, "Expecting the `Species` section in your YAML or JSON database to be either a list or dictionary. Please check other Reaktoro databases in either YAML or JSON format and replicate the structure.");

        }
    }

    /// Return the Data object with the details of an element with given unique @p symbol.
    auto getElementDetails(String const& symbol) -> Data
    {
        if(doc.exists("Elements"))
            if(doc["Elements"].exists(symbol))
                return doc["Elements"][symbol];
        return {};
    }

    /// Return the Data object with the details of a species with given unique @p name.
    auto getSpeciesDetails(const String& name) -> Data
    {
        if(doc.exists("Species"))
            if(doc["Species"].exists(name))
                return doc["Species"][name];
        return {};
    }

    /// Add a new element with given symbol. Check first if a Data object with key equals to element symbol exists.
    auto addElement(String const& symbol) -> Element
    {
        const auto attributes = getElementDetails(symbol);
        if(!attributes.isNull())
            return addElement(symbol, attributes); // create an element using info in the database file.
        Element element(symbol); // create an element using info from default elements in Elements.
        element_list.append(element);
        return element;
    }

    /// Add a new element with given `symbol` and `attributes`.
    auto addElement(String const& symbol, Data const& attributes) -> Element
    {
        errorif(!attributes.isDict(), "Expecting the attributes of an element as a dictionary in the Data object, but got instead:\n\n", attributes.repr());
        Element::Attribs attribs;
        attribs.symbol = symbol;
        attributes.at("MolarMass").to(attribs.molar_mass);
        if(attributes.exists("Name")) attributes.at("Name").to(attribs.name);
        if(attributes.exists("Tags")) attributes.at("Tags").to(attribs.tags);
        Element element(attribs);
        element_list.append(element);
        return element;
    }

    /// Add a new species with given unique @p name. A Data object for this species must exist.
    auto addSpecies(String const& name) -> Species
    {
        const auto attributes = getSpeciesDetails(name);
        errorif(attributes.isNull(), "Could not create a Species object with "
            "name `", name, "`, which does not seem to exist in the database. "
            "Are you sure this name is correct and there is a species with this name in the database?");
        return addSpecies(name, attributes);
    }

    /// Add a new species with given `name` and `attributes`.
    auto addSpecies(String const& name, Data const& attributes) -> Species
    {
        errorif(!attributes.isDict(), "Expecting the attributes of a species as an object, but got instead:\n\n", attributes.repr());
        errorif(!attributes.exists("Formula"), "Missing `Formula` specification in:\n\n", attributes.repr());
        errorif(!attributes.exists("AggregateState"), "Missing `AggregateState` specification in:\n\n", attributes.repr());
        errorif(!attributes.exists("Elements"), "Missing `Elements` specification in:\n\n", attributes.repr(), "\n",
            "Please assign `Elements: null` if this species does not have chemical elements (e.g., e-, which may be represented with only `Charge: -1`).");
        errorif(!attributes.exists("FormationReaction") && !attributes.exists("StandardThermoModel"), "Missing `FormationReaction` or `StandardThermoModel` specification in:\n\n", attributes.repr());
        const auto idx = species_list.find(name);
        if(idx < species_list.size())
            return species_list[idx]; // Do not add a species that has already been added! Return existing one.
        Species::Attribs attribs;
        attribs.name = name;
        attributes.at("Formula").to(attribs.formula);
        if(attributes.exists("Substance")) attributes.at("Substance").to(attribs.substance);
        if(attributes.exists("Charge")) attributes.at("Charge").to(attribs.charge);
        attributes.at("AggregateState").to(attribs.aggregate_state);
        errorif(attribs.aggregate_state == AggregateState::Undefined,
            "Unsupported AggregateState value `", attributes["AggregateState"].asString(), "` in:\n\n", attributes.repr(), "\n\n"
            "The supported values are given below:\n\n", supportedAggregateStateValues());
        attribs.elements = createElementalComposition(attributes);
        attribs.formation_reaction = createFormationReaction(attributes);
        attribs.std_thermo_model = createStandardThermoModel(attributes);
        attribs.tags = createTags(attributes);
        Species species(attribs);
        species_list.append(species);
        return species;
    }

    /// Create the elemental composition of the species whose attributes are found in the given Data object `data`.
    auto createElementalComposition(Data const& attributes) -> ElementalComposition
    {
        if(!attributes.exists("Elements") || attributes.at("Elements").isNull())
            return {}; // for example, e- may be described with only Charge: -1 and Elements: null

        Pairs<Element, double> pairs;
        const auto elements = attributes["Elements"].asString();
        const auto symbols_and_coeffs = parseNumberStringPairs(elements);
        for(const auto& [symbol, coeff] : symbols_and_coeffs)
        {
            const auto idx = element_list.find(symbol);
            if(idx < element_list.size())
                pairs.emplace_back(element_list[idx], coeff);
            else
            {
                const auto new_element = addElement(symbol);
                pairs.emplace_back(new_element, coeff);
            }
        }
        return ElementalComposition(pairs);
    }

    /// Create the vector of tags with given Data object with attributes for a species.
    auto createTags(Data const& attributes) -> Strings
    {
        if(attributes.exists("Tags"))
        {
            if(attributes["Tags"].isString())
                return split(attributes["Tags"].asString());
            if(attributes["Tags"].isList())
                return attributes["Tags"].as<Strings>();
        }
        return {};
    }

    /// Create a standard thermodynamic model with given Data object with attributes for a species.
    auto createStandardThermoModel(Data const& attributes) -> StandardThermoModel
    {
        if(attributes.exists("StandardThermoModel"))
            return attributes.at("StandardThermoModel").as<StandardThermoModel>();
        return {};
    }

    /// Create a formation reaction with given given Data object with attributes for a species.
    auto createFormationReaction(Data const& attributes) -> FormationReaction
    {
        if(attributes.exists("FormationReaction"))
            return FormationReaction()
                .withReactants(createReactants(attributes.at("FormationReaction")))
                .withReactionStandardThermoModel(createReactionStandardThermoModel(attributes.at("FormationReaction")))
                // .withProductStandardVolumeModel(createStandardVolumeModel(attributes.at("FormationReaction")))
                ;
        return {};
    }

    /// Create the list of reactant species in given formation reaction as a Data object.
    auto createReactants(Data const& data) -> Pairs<Species, double>
    {
        errorif(!data.exists("Reactants"), "Missing `Reactants` specification in:\n\n", data.repr());
        const auto names_and_coeffs = parseNumberStringPairs(data["Reactants"].asString());
        Pairs<Species, double> reactants;
        for(const auto& [name, coeff] : names_and_coeffs)
        {
            const auto idx = species_list.find(name);
            if(idx < species_list.size()) // check if there is a Species object with current name
                reactants.emplace_back(species_list[idx], coeff); // if so, add it to the list of reactants
            else // otherwise, add a new Species object
            {
                const auto new_species = addSpecies(name); // Note potential recursivity: addSpecies may also need to call createReactants! This is needed in case reactions are defined recursively.
                reactants.emplace_back(new_species, coeff);
            }
        }
        return reactants;
    }

    /// Create a standard thermodynamic model with given formation reaction as a Data object.
    auto createReactionStandardThermoModel(Data const& data) -> ReactionStandardThermoModel
    {
        errorif(!data.exists("ReactionStandardThermoModel"), "Missing `ReactionStandardThermoModel` specification in:\n\n", data.repr());
        return data["ReactionStandardThermoModel"].as<ReactionStandardThermoModel>();
    }

    // auto createStandardVolumeModel(Data const& child) -> StandardVolumeModel
    // {

    // }
};

DatabaseParser::DatabaseParser()
: pimpl(new Impl())
{}

DatabaseParser::DatabaseParser(const DatabaseParser& other)
: pimpl(new Impl(*other.pimpl))
{}

DatabaseParser::DatabaseParser(Data const& doc)
: pimpl(new Impl(doc))
{}

DatabaseParser::DatabaseParser(Data const& doc, const ElementList& elements)
: pimpl(new Impl(doc, elements))
{}

DatabaseParser::~DatabaseParser()
{}

auto DatabaseParser::operator=(DatabaseParser other) -> DatabaseParser&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto DatabaseParser::elements() const -> const ElementList&
{
    return pimpl->element_list;
}

auto DatabaseParser::species() const -> const SpeciesList&
{
    return pimpl->species_list;
}

DatabaseParser::operator Database() const
{
    return Database(elements(), species());
}

} // namespace Reaktoro
