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
using namespace Catch;

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/YAML.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Models/ReactionThermoModelConstLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelGemsLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelPhreeqcLgK.hpp>
#include <Reaktoro/Models/ReactionThermoModelVantHoff.hpp>
#include <Reaktoro/Models/StandardThermoModelConstant.hpp>
#include <Reaktoro/Models/StandardThermoModelHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelHollandPowell.hpp>
#include <Reaktoro/Models/StandardThermoModelInterpolation.hpp>
#include <Reaktoro/Models/StandardThermoModelMaierKelley.hpp>
#include <Reaktoro/Models/StandardThermoModelMineralHKF.hpp>
#include <Reaktoro/Models/StandardThermoModelWaterHKF.hpp>
#include <Reaktoro/Models/StandardVolumeModelConstant.hpp>
#include <Reaktoro/Serialization/Models.YAML.hpp>
using namespace Reaktoro;

//======================================================================
// Parameters in YAML format for standard thermodynamic models (stm)
//======================================================================
String params_stm_const = R"(
G0: 1.0
H0: 2.0
V0: 3.0
Cp0: 4.0
Cv0: 5.0
)";

String params_stm_mk = R"(
Gf: -3679250.6
Hf: -3876463.4
Sr: 209.32552
Vr: 9.281e-05
a: 251.41656
b: 0.0476976
c: -4769760.0
Tmax: 1700.0
)";

String params_stm_hkfmk1 = R"(
Gf: -3708312.7
Hf: -3931621.1
Sr: 207.14984
Vr: 0.00010025
ntr: 1
a: [258.1528, 342.58592]
b: [0.0581576, 0.014869936]
c: [-6280184.0, -20984434.0]
Ttr: [473.0]
Htr: [.nan]
Vtr: [.nan]
dPdTtr: [.nan]
Tmax: 1200.0
)";

String params_stm_hkfmk2 = R"(
Gf: -39522.064
Hf: -31589.2
Sr: 143.5112
Vr: 3.42e-05
ntr: 2
a: [65.39592, 7.610696, 90.3744]
b: [0.0359824, 0.221752, 0.0]
c: [0.0, 0.0, 0.0]
Ttr: [450.0, 620.0]
Htr: [3974.8, 2510.4]
Vtr: [.nan, .nan]
dPdTtr: [.nan, .nan]
Tmax: 1000.0
)";

String params_stm_hkfmk3 = R"(
Gf: .nan
Hf: .nan
Sr: 286.604
Vr: 0.0001432
ntr: 3
a: [369.61456, 409.65544, 488.18912, 461.24416]
b: [0.22643808, 0.14786256, 0.03112896, 0.04217472]
c: [-7556304.0, -4167264.0, -1937192.0, -1937192.0]
Ttr: [848.0, 950.0, 1050.0]
Htr: [.nan, .nan, .nan]
Vtr: [.nan, .nan, .nan]
dPdTtr: [.nan, .nan, .nan]
Tmax: 1100.0
)";

String params_stm_whkf = R"(
Ttr: 273.16
Str: 63.312288
Gtr: -235517.36
Htr: -287721.128
)";

String params_stm_hkf = R"(
Gf: 39371.44
Hf: -151084.24
Sr: 197.4848
a1: 5.8268894e-05
a2: 8251.0572
a3: 0.00049988758
a4: -150377.14
c1: 384.55521
c2: 116047.42
wref: -156816.32
charge: 0.0
Tmax: 0.0
)";

String params_stm_hp = R"(
Gf: -4937500.0
Hf: -5260650.0
Sr: 342.0
Vr: 0.00011525
a: 677.3
b: 0.0
c: -3772700.0
d: -5044.0
alpha0: 2.12e-05
kappa0: 190000000000.0
kappa0p: 2.98
kappa0pp: -1.6e-11
numatoms: 20.0
Tmax: 0.0
)";

String params_stm_hpg = R"(
Gf: -50710.0
Hf: -74810.0
Sr: 186.26
Vr: 0.0
a: 150.1
b: 0.002063
c: 3427700.0
d: -2650.4
Tmax: 0.0
)";

String params_stm_hpl = R"(
Gf: -2192340.0
Hf: -2307040.0
Sr: 127.6
Vr: 5.16e-05
a: 247.5
b: -0.003206
c: 0.0
d: -2051.9
alpha0: 2.9e-05
kappa0: 98500000000.0
kappa0p: 4.07
kappa0pp: -4.1e-11
numatoms: 7.0
Tmax: 0.0
)";

String params_stm_interpolation = R"(
Temperatures: [100, 200, 300]
Pressures: [400, 500]
G0: [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
H0: [[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
)";

//======================================================================
// Parameters in YAML format for standard volume models (svm)
//======================================================================
String params_svm_const = R"(
V0: 1.23e-5
)";

template<typename VecType>
auto allNaN(const VecType& v)
{
    for(auto&& x : v)
        if(!std::isnan(static_cast<double>(x))) return false;
    return true;
}

TEST_CASE("Testing Serialization for StandardThermoModel types", "[Serialization]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();

    SECTION("Testing YAML serialization of StandardThermoModelParamsConstant")
    {
        StandardThermoModelParamsConstant params = yaml::parse(params_stm_const);
        CHECK( params.G0  == 1.0 );
        CHECK( params.H0  == 2.0 );
        CHECK( params.V0  == 3.0 );
        CHECK( params.Cp0 == 4.0 );
        CHECK( params.Cv0 == 5.0 );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsMaierKelley")
    {
        StandardThermoModelParamsMaierKelley params = yaml::parse(params_stm_mk);
        CHECK( params.Gf   == -3679250.6 );
        CHECK( params.Hf   == -3876463.4 );
        CHECK( params.Sr   ==  209.32552 );
        CHECK( params.Vr   ==  9.281e-05 );
        CHECK( params.a    ==  251.41656 );
        CHECK( params.b    ==  0.0476976 );
        CHECK( params.c    == -4769760.0 );
        CHECK( params.Tmax ==  1700.0    );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsMineralHKF")
    {
        StandardThermoModelParamsMineralHKF params = yaml::parse(params_stm_hkfmk1);
        CHECK( params.Gf     == -3708312.7   );
        CHECK( params.Hf     == -3931621.1   );
        CHECK( params.Sr     ==  207.14984   );
        CHECK( params.Vr     ==  0.00010025  );
        CHECK( params.a[0]   ==  258.1528    );
        CHECK( params.a[1]   ==  342.58592   );
        CHECK( params.b[0]   ==  0.0581576   );
        CHECK( params.b[1]   ==  0.014869936 );
        CHECK( params.c[0]   == -6280184.0   );
        CHECK( params.c[1]   == -20984434.0  );
        CHECK( params.Ttr[0] ==  473.0       );
        CHECK( allNaN( params.Htr )          );
        CHECK( allNaN( params.Vtr )          );
        CHECK( allNaN( params.dPdTtr )       );
        CHECK( params.Tmax   == 1200.0       );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsMineralHKF")
    {
        StandardThermoModelParamsMineralHKF params = yaml::parse(params_stm_hkfmk2);
        CHECK( params.Gf     == -39522.064 );
        CHECK( params.Hf     == -31589.2   );
        CHECK( params.Sr     ==  143.5112  );
        CHECK( params.Vr     ==  3.42e-05  );
        CHECK( params.a[0]   ==  65.39592  );
        CHECK( params.a[1]   ==  7.610696  );
        CHECK( params.a[2]   ==  90.3744   );
        CHECK( params.b[0]   ==  0.0359824 );
        CHECK( params.b[1]   ==  0.221752  );
        CHECK( params.b[2]   ==  0.0       );
        CHECK( params.c[0]   ==  0.0       );
        CHECK( params.c[1]   ==  0.0       );
        CHECK( params.c[2]   ==  0.0       );
        CHECK( params.Ttr[0] ==  450.0     );
        CHECK( params.Ttr[1] ==  620.0     );
        CHECK( params.Htr[0] ==  3974.8    );
        CHECK( params.Htr[1] ==  2510.4    );
        CHECK( allNaN( params.Vtr )        );
        CHECK( allNaN( params.dPdTtr )     );
        CHECK( params.Tmax   == 1000.0     );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsMineralHKF")
    {
        StandardThermoModelParamsMineralHKF params = yaml::parse(params_stm_hkfmk3);
        CHECK( std::isnan( params.Gf )     );
        CHECK( std::isnan( params.Hf )     );
        CHECK( params.Sr     == 286.604    );
        CHECK( params.Vr     == 0.0001432  );
        CHECK( params.a[0]   == 369.61456  );
        CHECK( params.a[1]   == 409.65544  );
        CHECK( params.a[2]   == 488.18912  );
        CHECK( params.a[3]   == 461.24416  );
        CHECK( params.b[0]   == 0.22643808 );
        CHECK( params.b[1]   == 0.14786256 );
        CHECK( params.b[2]   == 0.03112896 );
        CHECK( params.b[3]   == 0.04217472 );
        CHECK( params.c[0]   == -7556304.0 );
        CHECK( params.c[1]   == -4167264.0 );
        CHECK( params.c[2]   == -1937192.0 );
        CHECK( params.c[3]   == -1937192.0 );
        CHECK( params.Ttr[0] == 848.0      );
        CHECK( params.Ttr[1] == 950.0      );
        CHECK( params.Ttr[2] == 1050.0     );
        CHECK( allNaN( params.Htr )        );
        CHECK( allNaN( params.Vtr )        );
        CHECK( allNaN( params.dPdTtr )     );
        CHECK( params.Tmax == 1100.0       );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsWaterHKF")
    {
        StandardThermoModelParamsWaterHKF params = yaml::parse(params_stm_whkf);
        CHECK( params.Ttr == 273.16 );
        CHECK( params.Str == 63.312288 );
        CHECK( params.Gtr == -235517.36 );
        CHECK( params.Htr == -287721.128 );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsHKF")
    {
        StandardThermoModelParamsHKF params = yaml::parse(params_stm_hkf);
        CHECK( params.Gf     == 39371.44      );
        CHECK( params.Hf     == -151084.24    );
        CHECK( params.Sr     == 197.4848      );
        CHECK( params.a1     == 5.8268894e-05 );
        CHECK( params.a2     == 8251.0572     );
        CHECK( params.a3     == 0.00049988758 );
        CHECK( params.a4     == -150377.14    );
        CHECK( params.c1     == 384.55521     );
        CHECK( params.c2     == 116047.42     );
        CHECK( params.wref   == -156816.32    );
        CHECK( params.charge == 0.0           );
        CHECK( params.Tmax   == 0.0           );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsHollandPowell")
    {
        StandardThermoModelParamsHollandPowell params = yaml::parse(params_stm_hp);
        CHECK( params.Gf       == -4937500.0     );
        CHECK( params.Hf       == -5260650.0     );
        CHECK( params.Sr       == 342.0          );
        CHECK( params.Vr       == 0.00011525     );
        CHECK( params.a        == 677.3          );
        CHECK( params.b        == 0.0            );
        CHECK( params.c        == -3772700.0     );
        CHECK( params.d        == -5044.0        );
        CHECK( params.alpha0   == 2.12e-05       );
        CHECK( params.kappa0   == 190000000000.0 );
        CHECK( params.kappa0p  == 2.98           );
        CHECK( params.kappa0pp == -1.6e-11       );
        CHECK( params.numatoms == 20.0           );
        CHECK( params.Tmax     == 0.0            );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsHollandPowell")
    {
        StandardThermoModelParamsHollandPowell params = yaml::parse(params_stm_hpg);
        CHECK( params.Gf       == -50710.0  );
        CHECK( params.Hf       == -74810.0  );
        CHECK( params.Sr       == 186.26    );
        CHECK( params.Vr       == 0.0       );
        CHECK( params.a        == 150.1     );
        CHECK( params.b        == 0.002063  );
        CHECK( params.c        == 3427700.0 );
        CHECK( params.d        == -2650.4   );
        CHECK( params.alpha0   == 0.0       );
        CHECK( params.kappa0   == 0.0       );
        CHECK( params.kappa0p  == 0.0       );
        CHECK( params.kappa0pp == 0.0       );
        CHECK( params.numatoms == 0.0       );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsHollandPowell")
    {
        StandardThermoModelParamsHollandPowell params = yaml::parse(params_stm_hpl);
        CHECK( params.Gf       == -2192340.0    );
        CHECK( params.Hf       == -2307040.0    );
        CHECK( params.Sr       == 127.6         );
        CHECK( params.Vr       == 5.16e-05      );
        CHECK( params.a        == 247.5         );
        CHECK( params.b        == -0.003206     );
        CHECK( params.c        == 0.0           );
        CHECK( params.d        == -2051.9       );
        CHECK( params.alpha0   == 2.9e-05       );
        CHECK( params.kappa0   == 98500000000.0 );
        CHECK( params.kappa0p  == 4.07          );
        CHECK( params.kappa0pp == -4.1e-11      );
        CHECK( params.numatoms == 7.0           );
        CHECK( params.Tmax     == 0.0           );
        // CHECK( params.Tcr      == 1710.0        );
        // CHECK( params.Smax     == 10.03         );
        // CHECK( params.Vmax     == 5.0e-07       );
    }

    SECTION("Testing YAML serialization of StandardThermoModelParamsInterpolation")
    {
        StandardThermoModelParamsInterpolation params = yaml::parse(params_stm_interpolation);
        CHECK( params.temperatures == Vec<double>{100, 200, 300} );
        CHECK( params.pressures == Vec<double>{400, 500} );
        CHECK( params.G0 == Vec<Vec<double>>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}} );
        CHECK( params.H0 == Vec<Vec<double>>{{4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}} );
        CHECK( params.V0.empty() );
        CHECK( params.Cp0.empty() );
        CHECK( params.Cv0.empty() );
        CHECK( params.VT0.empty() );
        CHECK( params.VP0.empty() );
    }
}

TEST_CASE("Testing Serialization for StandardVolumeModel types", "[Serialization]")
{
    SECTION("Testing YAML serialization of StandardVolumeModelParamsConstant")
    {
        StandardVolumeModelParamsConstant params = yaml::parse(params_svm_const);

        CHECK( params.V0 == 1.23e-5 );
    }
}
