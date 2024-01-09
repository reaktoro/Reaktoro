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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelMaierKelley.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StandardThermoModelMaierKelley class", "[StandardThermoModelMaierKelley]")
{
    const auto T = 75.0 + 273.15; // 75 degC (in K)
    const auto P = 1000.0 * 1e5;  // 1kbar (in Pa)

    // Check Oelkers et al. (1995), page 1553, table for CO2(g).
    SECTION("testing standard thermodynamic properties for CO2(g)")
    {
        // Parameters for CO2(g) from slop98.dat (converted to SI units)
        StandardThermoModelParamsMaierKelley params;
        params.Gf   = -394358.74;
        params.Hf   = -393509.38;
        params.Sr   =  213.73964;
        params.Vr   =  0.0;
        params.a    =  44.22488;
        params.b    =  0.0087864;
        params.c    = -861904.0;
        params.Tmax =  2500.0;

        auto model = StandardThermoModelMaierKelley(params);

        //======================================================================
        // Test method Model::operator()(T, P)
        //======================================================================

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-405198.0) ); // converted to J/mol from -96.84 kcal/mol as in table
        CHECK( props.H0  == Approx(-391571.0) );
        CHECK( props.V0  == Approx(0.0)       );
        CHECK( props.VT0 == Approx(0.0)       );
        CHECK( props.VP0 == Approx(0.0)       );
        CHECK( props.Cp0 == Approx(40.1729)   );

        //======================================================================
        // Test method Model::params()
        //======================================================================

        CHECK( model.params().isDict() );
        CHECK( model.params().at("MaierKelley").at("Gf").asFloat() == params.Gf );
        CHECK( model.params().at("MaierKelley").at("Hf").asFloat() == params.Hf );
        CHECK( model.params().at("MaierKelley").at("Sr").asFloat() == params.Sr );
        CHECK( model.params().at("MaierKelley").at("Vr").asFloat() == params.Vr );
        CHECK( model.params().at("MaierKelley").at("a").asFloat()  == params.a );
        CHECK( model.params().at("MaierKelley").at("b").asFloat()  == params.b );
        CHECK( model.params().at("MaierKelley").at("c").asFloat()  == params.c );
    }

    // Check Oelkers et al. (1995), page 1553, table for NH3(g).
    SECTION("testing standard thermodynamic properties for NH3(g)")
    {
        // Parameters for NH3(g) from slop98.dat (converted to SI units)
        StandardThermoModelParamsMaierKelley params;
        params.Gf   = -16451.488;
        params.Hf   = -46111.864;
        params.Sr   =  192.464;
        params.Vr   =  0.0;
        params.a    =  29.74824;
        params.b    =  0.025104;
        params.c    = -154808.0;
        params.Tmax =  1800.0;

        auto model = StandardThermoModelMaierKelley(params);

        StandardThermoProps props;
        props = model(T, P);

        CHECK( props.G0  == Approx(-26218.1) ); // converted to J/mol from -6.27 kcal/mol as in table
        CHECK( props.H0  == Approx(-44293.4) );
        CHECK( props.V0  == Approx(0.0)      );
        CHECK( props.VT0 == Approx(0.0)      );
        CHECK( props.VP0 == Approx(0.0)      );
        CHECK( props.Cp0 == Approx(37.211)   );
    }
}
