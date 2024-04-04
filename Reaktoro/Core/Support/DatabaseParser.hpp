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
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

// Forward declarations
class Database;
class Data;

/// Used to handle the parsing of YAML or JSON files to construct a Database object.
class DatabaseParser
{
public:
    /// Construct a default DatabaseParser object.
    DatabaseParser();

    /// Construct a copy of a DatabaseParser object.
    DatabaseParser(const DatabaseParser& other);

    /// Construct a DatabaseParser object with given Data object.
    explicit DatabaseParser(const Data& node);

    /// Construct a DatabaseParser object with given Data and Element objects.
    /// Use this if there are Element objects that should be reused to create the Species objects.
    explicit DatabaseParser(const Data& node, const ElementList& elements);

    /// Destroy this DatabaseParser object.
    ~DatabaseParser();

    /// Assign another DatabaseParser object to this.
    auto operator=(DatabaseParser other) -> DatabaseParser&;

    /// Return the parsed Element objects in the database file.
    auto elements() const -> const ElementList&;

    /// Return the parsed Species objects in the database file.
    auto species() const -> const SpeciesList&;

    /// Return the parsed Element objects in the database file.
    operator Database() const;

private:
    struct Impl;

    Ptr<Impl> pimpl;
};

} // namespace Reaktoro
