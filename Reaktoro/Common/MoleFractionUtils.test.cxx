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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/MoleFractionUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing MoleFractionUtils module", "[MoleFractionUtils]")
{
    // The possibilities for an array of species amounts
    auto n           = ArrayXr{{1.0, 2.0, 3.0, 4.0}};
    auto nzero       = ArrayXr{{0.0, 0.0, 0.0, 0.0}};
    auto nsingle     = ArrayXr{{2.0}};
    auto nsinglezero = ArrayXr{{0.0}};

    ArrayXr x, lnx;
    MatrixXd dxdn, dlnxdn;

    //-------------------------------------------------------------------------
    // TESTING METHOD: moleFractions
    //-------------------------------------------------------------------------

    x = moleFractions(n);

    CHECK( x[0] == Approx(0.1) );
    CHECK( x[1] == Approx(0.2) );
    CHECK( x[2] == Approx(0.3) );
    CHECK( x[3] == Approx(0.4) );

    x = moleFractions(nzero);

    CHECK( x[0] == Approx(0.0) );
    CHECK( x[1] == Approx(0.0) );
    CHECK( x[2] == Approx(0.0) );
    CHECK( x[3] == Approx(0.0) );

    x = moleFractions(nsingle);

    CHECK( x[0] == Approx(1.0) );

    x = moleFractions(nsinglezero);

    CHECK( x[0] == Approx(1.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: moleFractionsJacobian
    //-------------------------------------------------------------------------
    dxdn = moleFractionsJacobian(n);

    const auto dxdn_expected = MatrixXd{{
        { 0.09, -0.01, -0.01, -0.01},
        {-0.02,  0.08, -0.02, -0.02},
        {-0.03, -0.03,  0.07, -0.03},
        {-0.04, -0.04, -0.04,  0.06}
    }};

    INFO("dxdn = \n" << dxdn);
    CHECK( dxdn.isApprox(MatrixXd{{
        { 0.09, -0.01, -0.01, -0.01},
        {-0.02,  0.08, -0.02, -0.02},
        {-0.03, -0.03,  0.07, -0.03},
        {-0.04, -0.04, -0.04,  0.06},
    }}));

    dxdn = moleFractionsJacobian(nzero);

    INFO("dxdn = \n" << dxdn);
    CHECK( dxdn.isApprox(MatrixXd::Zero(4, 4)) );

    dxdn = moleFractionsJacobian(nsingle);

    INFO("dxdn = \n" << dxdn);
    CHECK( dxdn.isApprox(MatrixXd::Zero(1, 1)) );

    dxdn = moleFractionsJacobian(nsinglezero);

    INFO("dxdn = \n" << dxdn);
    CHECK( dxdn.isApprox(MatrixXd::Zero(1, 1)) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: lnMoleFractionsJacobian
    //-------------------------------------------------------------------------
    dlnxdn = lnMoleFractionsJacobian(n);

    INFO("dlnxdn = \n" << dlnxdn);
    CHECK( dlnxdn.isApprox(MatrixXd{{
        { 0.9/1,   -0.1,   -0.1,   -0.1},
        {  -0.1,  0.8/2,   -0.1,   -0.1},
        {  -0.1,   -0.1,  0.7/3,   -0.1},
        {  -0.1,   -0.1,   -0.1,  0.6/4},
    }}));

    dlnxdn = lnMoleFractionsJacobian(nzero);

    INFO("dlnxdn = \n" << dlnxdn);
    CHECK( dlnxdn.isApprox(MatrixXd::Zero(4, 4)) );

    dlnxdn = lnMoleFractionsJacobian(nsingle);

    INFO("dlnxdn = \n" << dlnxdn);
    CHECK( dlnxdn.isApprox(MatrixXd::Zero(1, 1)) );

    dlnxdn = lnMoleFractionsJacobian(nsinglezero);

    INFO("dlnxdn = \n" << dlnxdn);
    CHECK( dlnxdn.isApprox(MatrixXd::Zero(1, 1)) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: lnMoleFractionsJacobianDiagonal
    //-------------------------------------------------------------------------
    ArrayXd diag(4);
    ArrayXd diag_expected(4);

    diag = lnMoleFractionsJacobianDiagonal(n);

    diag_expected = 1/n - 1/n.sum();

    INFO("diag = \n" << diag);
    INFO("diag(expected) = \n" << diag_expected);
    CHECK( diag.isApprox(diag_expected) );
}
