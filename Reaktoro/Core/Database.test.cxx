// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/Database.hpp>
using namespace Reaktoro;

namespace test {

/// Return mock standard thermodynamic model for an aqueous species.
auto standardThermoModelAqueous(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 2.1 * log(P/T);
    props.H0  = 2.2 * log(P/T);
    props.V0  = 2.3 * log(P/T);
    props.Cp0 = 2.4 * log(P/T);
    return props;
}

/// Return mock standard thermodynamic model for a gaseous species.
auto standardThermoModelGaseous(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 0.1 * log(P/T);
    props.H0  = 0.2 * log(P/T);
    props.V0  = 0.3 * log(P/T);
    props.Cp0 = 0.4 * log(P/T);
    return props;
}

/// Return mock standard thermodynamic model for a solid species.
auto standardThermoModelSolid(real T, real P) -> StandardThermoProps
{
    StandardThermoProps props;
    props.G0  = 1.1 * log(P/T);
    props.H0  = 1.2 * log(P/T);
    props.V0  = 1.3 * log(P/T);
    props.Cp0 = 1.4 * log(P/T);
    return props;
}

/// Return a mock Species object of aqueous type for test reasons.
auto createAqueousSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Aqueous)
        .withStandardThermoModel(standardThermoModelAqueous);
}

/// Return a mock Species object of gas type for test reasons.
auto createGaseousSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Gas)
        .withStandardThermoModel(standardThermoModelGaseous);
}

/// Return a mock Species object of solid type for test reasons.
auto createSolidSpecies(String name)
{
    return Species(name)
        .withAggregateState(AggregateState::Solid)
        .withStandardThermoModel(standardThermoModelSolid);
}

/// Return a mock Database object for test reasons.
auto createDatabase() -> Database
{
    Database db;

    // AQUEOUS SPECIES
    db.addSpecies( createAqueousSpecies("H2O(aq)"  ) );
    db.addSpecies( createAqueousSpecies("H+(aq)"   ) );
    db.addSpecies( createAqueousSpecies("OH-(aq)"  ) );
    db.addSpecies( createAqueousSpecies("H2(aq)"   ) );
    db.addSpecies( createAqueousSpecies("O2(aq)"   ) );
    db.addSpecies( createAqueousSpecies("Na+(aq)"  ) );
    db.addSpecies( createAqueousSpecies("Cl-(aq)"  ) );
    db.addSpecies( createAqueousSpecies("NaCl(aq)" ) );
    db.addSpecies( createAqueousSpecies("HCl(aq)"  ) );
    db.addSpecies( createAqueousSpecies("NaOH(aq)" ) );
    db.addSpecies( createAqueousSpecies("Ca++(aq)" ) );
    db.addSpecies( createAqueousSpecies("Mg++(aq)" ) );
    db.addSpecies( createAqueousSpecies("CO2(aq)"  ) );
    db.addSpecies( createAqueousSpecies("HCO3-(aq)") );
    db.addSpecies( createAqueousSpecies("CO3--(aq)") );
    db.addSpecies( createAqueousSpecies("CaCl2(aq)") );
    db.addSpecies( createAqueousSpecies("MgCl2(aq)") );
    db.addSpecies( createAqueousSpecies("SiO2(aq)" ) );
    db.addSpecies( createAqueousSpecies("e-(aq)" ) );

    // GASEOUS SPECIES
    db.addSpecies( createGaseousSpecies("CO2(g)" ) );
    db.addSpecies( createGaseousSpecies("O2(g)"  ) );
    db.addSpecies( createGaseousSpecies("H2(g)"  ) );
    db.addSpecies( createGaseousSpecies("H2O(g)" ) );
    db.addSpecies( createGaseousSpecies("CH4(g)" ) );
    db.addSpecies( createGaseousSpecies("CO(g)"  ) );

    // SOLID SPECIES
    db.addSpecies( createSolidSpecies("NaCl(s)"      ) );
    db.addSpecies( createSolidSpecies("CaCO3(s)"     ) );
    db.addSpecies( createSolidSpecies("MgCO3(s)"     ) );
    db.addSpecies( createSolidSpecies("CaMg(CO3)2(s)") );
    db.addSpecies( createSolidSpecies("SiO2(s)"      ) );

    return db;
}

} // namespace test

TEST_CASE("Testing Database class", "[Database]")
{
    Database db = test::createDatabase();

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::elements
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.elements().index("H") );
    REQUIRE_NOTHROW( db.elements().index("C") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::species
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.species().index("H2O(aq)") );
    REQUIRE_NOTHROW( db.species().index("CO2(g)")  );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::speciesWithAggregateState
    //-------------------------------------------------------------------------
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H2O(aq)") );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Aqueous).index("H+(aq)")  );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Gas).index("CO2(g)")      );
    REQUIRE_NOTHROW( db.speciesWithAggregateState(AggregateState::Solid).index("NaCl(s)")   );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::addSpecies
    //-------------------------------------------------------------------------
    db.addSpecies( Species("Fe(s)").withName("Iron")   );
    db.addSpecies( Species("Cu(s)").withName("Copper") );

    REQUIRE_NOTHROW( db.species().index("Iron")   );
    REQUIRE_NOTHROW( db.species().index("Copper") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::attachData
    //-------------------------------------------------------------------------
    db.attachData( String("SomeData") );

    REQUIRE( std::any_cast<String>(db.attachedData()) == "SomeData" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::extends
    //-------------------------------------------------------------------------
    Database dbx;

    dbx.addSpecies( Species("H2O(aq)"  ) ); // check repeated species
    dbx.addSpecies( Species("CaCO3(aq)") );
    dbx.addSpecies( Species("MgCO3(aq)") );
    dbx.addSpecies( Species("Al+2(aq)" ) );
    dbx.addSpecies( Species("Al+3(aq)" ) );
    dbx.addSpecies( Species("Fe+2(aq)" ) );
    dbx.addSpecies( Species("Fe+3(aq)" ) );
    dbx.addSpecies( Species("Mg(s)") );
    dbx.addSpecies( Species("Al(s)") );

    db.extend(dbx);

    CHECK_NOTHROW( db.species().index("H2O(aq)!") ); // repeated species get ! added to name
    CHECK_NOTHROW( db.species().index("CaCO3(aq)") );
    CHECK_NOTHROW( db.species().index("MgCO3(aq)") );
    CHECK_NOTHROW( db.species().index("Al+2(aq)") );
    CHECK_NOTHROW( db.species().index("Al+3(aq)") );
    CHECK_NOTHROW( db.species().index("Fe+2(aq)") );
    CHECK_NOTHROW( db.species().index("Fe+3(aq)") );
    CHECK_NOTHROW( db.species().index("Mg(s)") );
    CHECK_NOTHROW( db.species().index("Al(s)") );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Database::clear
    //-------------------------------------------------------------------------
    db.clear();

    REQUIRE( db.elements().size() == 0 );
    REQUIRE( db.species().size()  == 0 );

    //-------------------------------------------------------------------------
    // TESTING CONSTRUCTORS: Database(species) and Database(elements, species)
    //-------------------------------------------------------------------------
    auto check_same_contents_in_databases = [](const Database& db1, const Database& db2)
    {
        CHECK( db1.species().size() == db2.species().size() );
        CHECK( db1.elements().size() == db2.elements().size() );

        for(auto obj : db1.elements())
            CHECK( containsfn(db2.elements(), RKT_LAMBDA(x, x.name() == obj.symbol())) );

        for(auto obj : db1.species())
            CHECK( containsfn(db2.species(), RKT_LAMBDA(x, x.name() == obj.name())) );
    };

    Database new_db1(db.species());
    Database new_db2(db.elements(), db.species());

    check_same_contents_in_databases(db, new_db1);
    check_same_contents_in_databases(db, new_db2);
}

TEST_CASE("Testing Database object creation using Database::fromContents", "[Database]")
{
    String contents = R"xyz(
Species:
- Name: Akermanite
  Formula: Ca2MgSi2O7
  Elements: 2:Ca 1:Mg 2:Si 7:O
  AggregateState: Solid
  StandardThermoModel:
    MaierKelley:
      Gf: -3679250.6
      Hf: -3876463.4
      Sr: 209.32552
      Vr: 9.281e-05
      a: 251.41656
      b: 0.0476976
      c: -4769760.0
      Tmax: 1700.0
    )xyz";

    Database db = Database::fromContents(contents);

    CHECK(db.species().size() == 1);
    CHECK(db.species()[0].name() == "Akermanite");
}
