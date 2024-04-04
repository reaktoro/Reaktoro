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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The class used to store and retrieve data of chemical species.
/// @see Element, Species
/// @ingroup Core
class Database
{
public:
    /// Return a Database object constructed with a given local file.
    /// @warning An exception is thrown if `path` does not point to a valid local database file.
    /// @param path The path, including file name, to the local database file.
    static auto fromFile(String const& path) -> Database;

    /// Return a Database object constructed with a given embedded file.
    /// @warning An exception is thrown if `path` does not point to a valid embedded database file.
    /// @param path The path, including file name, to the embedded database file.
    static auto fromEmbeddedFile(String const& path) -> Database;

    /// Return a Database object constructed with given database text contents.
    /// @param contents The contents of the database as a string.
    static auto fromContents(String const& contents) -> Database;

    /// Return a Database object constructed with given input stream containing the database text contents.
    /// @param stream The input stream containing the database file contents.
    static auto fromStream(std::istream& stream) -> Database;

    /// @copydoc Database::fromFile
    static auto local(String const& path) -> Database;

    /// @copydoc Database::fromEmbeddedFile
    static auto embedded(String const& path) -> Database;

    /// Construct a default Database object.
    Database();

    /// Construct a copy of a Database object.
    Database(Database const& other);

    /// Construct a Database object with given elements and species.
    Database(Vec<Element> const& elements, Vec<Species> const& species);

    /// Construct a Database object with given species (elements extracted from them).
    explicit Database(Vec<Species> const& species);

    /// Destroy this Database object.
    ~Database();

    /// Assign another Database object to this.
    auto operator=(Database other) -> Database&;

    /// Remove all species and elements from the database.
    auto clear() -> void;

    /// Add an element in the database.
    /// @note If an Element object with same symbol already exists in the
    /// Database container, the given Element object is not added.
    auto addElement(Element const& element) -> void;

    /// Add a species in the database.
    /// @note If a Species object with same name already exists in the Database
    /// container, the added Species object has its name slightly changed (e.g. `H2O` becomes `H2O!`).
    auto addSpecies(Species const& species) -> void;

    /// Add a list of species in the database.
    auto addSpecies(Vec<Species> const& species) -> void;

    /// Attach data to this database whose type is known at runtime only.
    auto attachData(Any const& data) -> void;

    /// Extend this database with elements, species and other contents from another database.
    auto extend(Database const& other) -> void;

    /// Extend this database with elements, species and other contents from another database.
    auto extendWithDatabase(Database const& other) -> void;

    /// Extend this database with elements, species and other contents from another database.
    /// @warning An exception is thrown if `path` does not point to a valid local database file.
    /// @param path The path, including file name, to the local database file.
    auto extendWithFile(String const& path) -> void;

    /// Return all elements in the database.
    auto elements() const -> ElementList const&;

    /// Return all species in the database.
    auto species() const -> SpeciesList const&;

    /// Return all species in the database with given aggregate state.
    auto speciesWithAggregateState(AggregateState option) const -> SpeciesList;

    /// Return an element with given symbol in the database.
    /// @warning An exception is thrown if no element with given symbol exists.
    auto element(String const& symbol) const -> Element const&;

    /// Return a species with given name in the database.
    /// @warning An exception is thrown if no species with given name exists.
    auto species(String const& name) const -> Species const&;

    /// Construct a reaction with given equation.
    /// @warning An exception is thrown if the reaction has an inexistent species in the database.
    auto reaction(String const& equation) const -> Reaction;

    /// Return the attached data to this database whose type is known at runtime only.
    auto attachedData() const -> Any const&;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
