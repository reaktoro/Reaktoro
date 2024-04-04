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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Element.hpp>

namespace Reaktoro {

/// A type used store a collection of elements.
/// @see Element
class Elements
{
public:
    /// Construct a copy of a Elements object [deleted].
    Elements(const Elements&) = delete;

    /// Assign a Elements object to this [deleted].
    auto operator=(const Elements&) -> Elements& = delete;

    /// Return the single Elements object.
    static auto instance() -> Elements&;

    /// Return the elements in the collection of elements.
    static auto data() -> Vec<Element> const&;

    /// Return the elements in the collection of elements.
    [[deprecated("Use Elements::data() instead.")]]
    static auto elements() -> Vec<Element> const&;

    /// Clear the elements in the collection of elements.
    static auto clear() -> void;

    /// Reset the elements in the collection of elements to its default state.
    static auto reset() -> void;

    /// Append a custom element to the collection of elements.
    static auto append(Element element) -> void;

    /// Replace an existing element in the collection of elements with another. If given existing element is not present in the collection, this method returns false.
    static auto replace(String const& existingSymbol, Element const& replacement) -> bool;

    /// Return the number of elements in the collection of elements.
    static auto size() -> std::size_t;

    /// Return the element with given symbol.
    static auto withSymbol(String symbol) -> Optional<Element>;

    /// Return the element with given name.
    static auto withName(String name) -> Optional<Element>;

    /// Return the elements with a given tag.
    static auto withTag(String tag) -> Vec<Element>;

    /// Return the elements with given tags.
    static auto withTags(StringList const& tags) -> Vec<Element>;

    /// Return begin const iterator of this Elements instance.
    auto begin() const { return data().begin(); }

    /// Return begin iterator of this Elements instance.
    auto begin() { return data().begin(); }

    /// Return end const iterator of this Elements instance.
    auto end() const { return data().end(); }

    /// Return end iterator of this Elements instance.
    auto end() { return data().end(); }

private:
    /// The elements stored in the collection of elements.
    Vec<Element> m_elements;

private:
    /// Construct a default Elements object [private].
    Elements();

    /// Destroy this Elements object [private].
    ~Elements();
};

} // namespace Reaktoro
