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

#include "Database.hpp"

// C++ includes
#include <fstream>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/ParseUtils.hpp>
#include <Reaktoro/Core/Data.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Core/Support/DatabaseParser.hpp>

namespace Reaktoro {

template<typename Source>
auto createDatabaseFromContents(Source& contents)
{
    Data doc;

    try { doc = Data::parseYaml(contents); }
    catch(...)
    {
        try { doc = Data::parseJson(contents); }
        catch(...)
        {
            errorif(true, "Could not parse given text in Database::fromContents as it does not seem to be either YAML or JSON formats.");
        }
    }

    DatabaseParser dbparser(doc);
    return Database(dbparser);
}

auto createDataFromJsonOrYaml(String const& path) -> Data
{
    std::ifstream file(path);
    errorif(!file.is_open(),
        "Could not open file `", path, "`. Ensure the given file path "
        "is relative to the directory where your application is RUNNING "
        "(not necessarily where the executable is located!). Alternatively, "
        "try a full path to the file (e.g., "
        "in Windows, `C:\\User\\username\\mydata\\mydatabase.yaml`, "
        "in Linux and macOS, `/home/username/mydata/mydatabase.yaml`). "
        "File formats accepted are JSON and YAML and expected file extensions are .json, .yaml, or .yml.");
    auto isJson = endswith(path, ".json");
    auto isYaml = endswith(path, ".yaml") || endswith(path, ".yml");
    errorifnot(isJson || isYaml, "The file `", path, "` must be a JSON or YAML file terminating with .json, .yaml, or .yml.");
    auto doc = isJson ? Data::parseJson(file) : Data::parseYaml(file);
    return doc;
}

struct Database::Impl
{
    /// The Species objects in the database.
    SpeciesList species;

    /// The Element objects in the database.
    ElementList elements;

    /// The additional data in the database whose type is known at runtime only.
    Any attached_data;

    /// The symbols of all elements already in the database.
    Set<String> inserted_elements_set;

    /// The ids of the species already inserted in the database, with id = species name + aggregate state
    Set<String> inserted_species_set;

    /// The species in the database grouped in terms of their aggregate state
    Map<AggregateState, SpeciesList> species_with_aggregate_state;

    /// Add an element in the database.
    auto addElement(Element const& element) -> void
    {
        auto const ielement = elements.findWithSymbol(element.symbol());

        if(ielement == elements.size())
        {
            elements.append(element);
            inserted_elements_set.insert(element.symbol());
        }
        else
        {
            errorif(element.molarMass() != elements[ielement].molarMass(), "Element with symbol `", element.symbol(), "` already exists in the database with molar mass ", elements[ielement].molarMass(), "kg/mol. You are trying to add a new element with molar mass ", element.molarMass(), "kg/mol. It's possible that your database has a duplicated element entry with inconsistent molar mass and other attributes");
            errorif(element.name() != elements[ielement].name(), "Element with symbol `", element.symbol(), "` already exists in the database with name `", elements[ielement].name(), "`. You are trying to add a new element with name `", element.name(), "`. It's possible that your database has a duplicated element entry with inconsistent names and other attributes");
        }
    }

    /// Add a species in the database.
    auto addSpecies(Species newspecies) -> void
    {
        // Ensure unique names for species with same aggregate state are used when storing a new species!
        auto const name = newspecies.name();
        auto const agstate = newspecies.aggregateState();

        // Find a unique name and id next, if needed.
        auto unique_id = name + std::to_string(static_cast<int>(agstate));
        auto unique_name = name;

        while(inserted_species_set.find(unique_id) != inserted_species_set.end())
        {
            unique_name = unique_name + "!"; // keep adding symbol ! to the name such as H2O, H2O!, H2O!!, H2O!!! if many H2O are given
            unique_id = unique_name + std::to_string(static_cast<int>(agstate));
        }

        // Replace original name with found unique name (using appended !) if name is not unique among the species with same aggregate state.
        if(name != unique_name) {
            newspecies = newspecies.withName(unique_name);
            warningif(true, "Species with same aggregate state should have unique names in Database, but species `", name, "` with aggregate state `", agstate, "` violates this rule. The unique name ", unique_name, " has been assigned instead.");
        }

        // Replace aggregate state to Aqueous if undefined. If this default
        // aggregate state option is not appropriate for a given scenario,
        // ensure then that the correct aggregate state of the species has
        // already been set in the Species object using method
        // Species::withAggregateState. This permits species names such as H2O,
        // CO2, H2, O2 as well as ions H+, OH-, HCO3- to be considered as
        // aqueous species, without requiring explicit aggregate state suffix
        // aq as in H2O(aq), CO2(aq), H2(aq). For all other species, ensure an
        // aggregate state suffix is provided, such as H2O(l), H2O(g),
        // CaCO3(s), Fe3O4(s).
        if(newspecies.aggregateState() == AggregateState::Undefined)
            newspecies = newspecies.withAggregateState(AggregateState::Aqueous);

        // Append the new Species object in the species container
        species.append(newspecies);

        // Update the list of unique species ids.
        inserted_species_set.insert(unique_id);

        // Update the container of elements with the Element objects in this new species.
        for(auto&& [element, coeff] : newspecies.elements())
            addElement(element);

        // Add the new species in the group of species with same aggregate state
        species_with_aggregate_state[newspecies.aggregateState()].push_back(newspecies);
    }

    /// Construct a reaction with given equation.
    auto reaction(String const& equation) const -> Reaction
    {
        return Reaction().withEquation(ReactionEquation(equation, species));
    }
};

Database::Database()
: pimpl(new Impl())
{}

Database::Database(Database const& other)
: pimpl(new Impl(*other.pimpl))
{}

Database::Database(Vec<Element> const& elements, Vec<Species> const& species)
: Database()
{
    for(auto const& x : elements)
        addElement(x);
    for(auto const& x : species)
        addSpecies(x);
}

Database::Database(Vec<Species> const& species)
: Database()
{
    for(auto const& x : species)
        addSpecies(x);
}

Database::~Database()
{}

auto Database::operator=(Database other) -> Database&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Database::clear() -> void
{
    *pimpl = Database::Impl();
}

auto Database::addElement(Element const& element) -> void
{
    pimpl->addElement(element);
}

auto Database::addSpecies(Species const& species) -> void
{
    pimpl->addSpecies(species);
}

auto Database::addSpecies(Vec<Species> const& species) -> void
{
    for(auto const& x : species)
        addSpecies(x);
}

auto Database::attachData(Any const& data) -> void
{
    pimpl->attached_data = data;
}

auto Database::extend(Database const& other) -> void
{
    extendWithDatabase(other);
}

auto Database::extendWithDatabase(Database const& other) -> void
{
    for(auto const& element : other.elements())
        addElement(element);

    for(auto const& species : other.species())
        addSpecies(species);

    // TODO: Replace Any by Map<String, Any> so that it becomes easier/more intuitive to unify different attached data to Database objects.
    // pimpl->attached_data = ???;
}

auto Database::extendWithFile(String const& path) -> void
{
    // Prvide current element objects in this database which can be reused to
    // create the Species objects in the extended database.
    auto doc = createDataFromJsonOrYaml(path);
    DatabaseParser dbparser(doc, pimpl->elements);
    Database dbx(dbparser);
    extendWithDatabase(dbx);
}

auto Database::elements() const -> ElementList const&
{
    return pimpl->elements;
}

auto Database::species() const -> SpeciesList const&
{
    return pimpl->species;
}

auto Database::speciesWithAggregateState(AggregateState option) const -> SpeciesList
{
    auto it = pimpl->species_with_aggregate_state.find(option);
    if(it == pimpl->species_with_aggregate_state.end())
        return {};
    return it->second;
}

auto Database::element(String const& symbol) const -> Element const&
{
    return elements().getWithSymbol(symbol);
}

auto Database::species(String const& name) const -> Species const&
{
    return species().getWithName(name);
}

auto Database::reaction(String const& equation) const -> Reaction
{
    return pimpl->reaction(equation);
}

auto Database::attachedData() const -> Any const&
{
    return pimpl->attached_data;
}

auto Database::fromFile(String const& path) -> Database
{
    auto doc = createDataFromJsonOrYaml(path);
    DatabaseParser dbparser(doc);
    return Database(dbparser);
}

auto Database::fromEmbeddedFile(String const& path) -> Database
{
    const String contents = Embedded::get("databases/reaktoro/" + path);
    return fromContents(contents);
}

auto Database::fromContents(String const& contents) -> Database
{
    return createDatabaseFromContents(contents);
}

auto Database::fromStream(std::istream& stream) -> Database
{
    return createDatabaseFromContents(stream);
}

auto Database::local(String const& path) -> Database
{
    return fromFile(path);
}

auto Database::embedded(String const& path) -> Database
{
    return fromEmbeddedFile(path);
}

} // namespace Reaktoro

