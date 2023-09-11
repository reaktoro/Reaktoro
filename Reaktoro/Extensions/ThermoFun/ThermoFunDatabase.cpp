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

#include "ThermoFunDatabase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Extensions/ThermoFun/ThermoFunEngine.hpp>

// ThermoFun includes
#include <ThermoFun/ThermoFun.h>

namespace Reaktoro {
namespace {

/// Return the standard thermodynamic property function of a species with given name.
auto createStandardThermoModel(const ThermoFunEngine& engine, const String& species) -> StandardThermoModel
{
    const auto& substances = engine.database().mapSubstances();
    const auto it = substances.find(species);

    errorif(it == substances.end(), "Expecting a species name that exists in the ThermoFun database, but got `", species, "` instead.");

    const auto substance = it->second;

    return [=](real T, real P) -> StandardThermoProps
    {
        return engine.props(T, P, substance);
    };
}

/// Convert a ThermoFun::Element object into a Reaktoro::Element object
auto createElement(const ThermoFun::Element& element) -> Element
{
    Element converted;
    converted = converted.withName(element.name());
    converted = converted.withSymbol(element.symbol());
    converted = converted.withMolarMass(element.molarMass() * 1e-3); // from g/mol to kg/mol
    return converted;
}

/// Return the elements and their coefficients in a species a ThermoFun::Substance object
auto createElements(const ThermoFunEngine& engine, const ThermoFun::Substance& substance) -> Pairs<Element, double>
{
    auto db = engine.database();
    Pairs<Element, double> elements;
    for(auto&& [element, coeff] : db.parseSubstanceFormula(substance.formula()))
    {
        if(element.symbol() == "Zz")
            continue; // skip charge element in ThermoFun (charge is explicitly accessed in Reaktoro::Species)
        elements.emplace_back(createElement(element), coeff);
    }
    return elements;
}

/// Convert aggregate state value from ThermoFun to Reaktoro::AggregateState
auto convertAggregateState(const ThermoFun::AggregateState::type& value)
{
    switch(value)
    {
    case ThermoFun::AggregateState::GAS:            return AggregateState::Gas;
    case ThermoFun::AggregateState::LIQUID:         return AggregateState::Liquid;
    case ThermoFun::AggregateState::GLASS:          return AggregateState::Vitreous;
    case ThermoFun::AggregateState::CRYSTAL:        return AggregateState::CrystallineSolid;
    case ThermoFun::AggregateState::AQUEOUS:        return AggregateState::Aqueous;
    case ThermoFun::AggregateState::SURFACE:        return AggregateState::Adsorbed;
    case ThermoFun::AggregateState::IONEX:          return AggregateState::IonExchange;
    case ThermoFun::AggregateState::PLASMA:         return AggregateState::Plasma;
    case ThermoFun::AggregateState::SOLID:          return AggregateState::Solid;
    case ThermoFun::AggregateState::CONDENSED:      return AggregateState::CondensedPhase;
    case ThermoFun::AggregateState::FLUID:          return AggregateState::Fluid;
    case ThermoFun::AggregateState::LIQUIDCRYSTAL:  return AggregateState::LiquidCrystal;
    case ThermoFun::AggregateState::AMORPHOUSSOLID: return AggregateState::AmorphousSolid;
    case ThermoFun::AggregateState::MONOMERIC:      return AggregateState::Monomeric;
    case ThermoFun::AggregateState::POLYMERIC:      return AggregateState::Polymeric;
    case ThermoFun::AggregateState::SOLIDSOLUTION:  return AggregateState::SolidSolution;
    case ThermoFun::AggregateState::OTHER:          return AggregateState::Undefined;
    default:                                        return AggregateState::Undefined;
    }
}

/// Convert a ThermoFun::Substance object into a Reaktoro::Species object
auto createSpecies(const ThermoFunEngine& engine, const ThermoFun::Substance& substance) -> Species
{
    Species species;
    species = species.withName(substance.symbol());
    species = species.withFormula(substance.formula());
    species = species.withSubstance(substance.name());
    species = species.withElements(createElements(engine, substance));
    species = species.withCharge(substance.charge());
    species = species.withAggregateState(convertAggregateState(substance.aggregateState()));
    species = species.withStandardThermoModel(createStandardThermoModel(engine, substance.symbol()));
    species = species.withAttachedData(substance);
    return species;
}

} // namespace

ThermoFunDatabase::ThermoFunDatabase()
{}

ThermoFunDatabase::ThermoFunDatabase(const String& name)
: ThermoFunDatabase(ThermoFunDatabase::withName(name))
{}

ThermoFunDatabase::ThermoFunDatabase(const ThermoFun::Database& db)
{
    ThermoFunEngine engine(db);
    attachData(engine);
    for(auto [_, subs] : db.mapSubstances())
        addSpecies(createSpecies(engine, subs));
}

auto ThermoFunDatabase::withName(const String& name) -> ThermoFunDatabase
{
    errorif(!oneof(name,
        "aq17",
        "cemdata18",
        "heracles",
        "mines16",
        "psinagra-12-07",
        "slop98-organic",
        "slop98"),
        "Could not load embedded ThermoFun database file with name `", name, "`. ",
        "The currently supported names are: \n"
        "    - aq17            (corresponding file: aq17-thermofun.json)      \n",
        "    - cemdata18       (corresponding file: cemdata18-thermofun.json) \n",
        "    - heracles        (corresponding file: heracles-thermofun.json)  \n",
        "    - mines16         (corresponding file: mines16-thermofun.json)   \n",
        "    - psinagra-12-07  (corresponding file: psinagra-thermofun.json)  \n",
        "    - slop98-organic  (corresponding file: slop98-thermofun.json)    \n",
        "    - slop98          (corresponding file: slop98-thermofun.json)    \n",
        "");
    const auto text = Embedded::get("databases/thermofun/" + name + "-thermofun.json");
    return fromFile(text);
}

auto ThermoFunDatabase::fromFile(const String& filepath) ->  ThermoFunDatabase
{
    ThermoFun::Database db(filepath);
    return ThermoFunDatabase(db);
}

auto ThermoFunDatabase::fromFiles(const Strings& filepaths) ->  ThermoFunDatabase
{
    errorif(filepaths.empty(), "Expecting at least one file path in method ThermoFunDatabase::fromFiles.");
    ThermoFun::Database db(filepaths[0]);
    for(auto i = 1; i < filepaths.size(); ++i)
        db.appendData(filepaths[i]);
    return ThermoFunDatabase(db);
}

auto ThermoFunDatabase::fromContents(const String& contents) ->  ThermoFunDatabase
{
    return ThermoFunDatabase::fromFile(contents);
}

} // namespace Reaktoro
