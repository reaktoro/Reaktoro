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
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Singletons/DissociationReactions.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
#include <Reaktoro/Water/WaterElectroProps.hpp>
#include <Reaktoro/Water/WaterElectroPropsJohnsonNorton.hpp>
#include <Reaktoro/Water/WaterThermoProps.hpp>
#include <Reaktoro/Water/WaterThermoPropsUtils.hpp>
using namespace Reaktoro;

auto moleFractions(Index size) -> ArrayXr
{
    const auto n = ArrayXr::Random(size);
    return n / n.sum();
}

TEST_CASE("Testing AqueousMixture class", "[AqueousMixture]")
{
    // Ensure the dissociation reactions are in its default state (which may have been changed by other tests!)
    DissociationReactions::reset();

    WHEN("When the aqueous mixture is setup correctly")
    {
        SpeciesList species("H2O H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO3-- K+ CO2 HCl NaCl NaOH CaCl2 MgCl2 CaCO3 MgCO3");

        AqueousMixture mixture(species);

        CHECK( mixture.species().size() == species.size() );

        // Test H2O is not classified as neutral solute species
        CHECK_THROWS( mixture.neutral().index("H2O") );

        for(auto x : species)
        {
            if(x.name() == "H2O") continue;

            // Test AqueousMixture::neutral|charged|cations|anions methods
            if(x.charge() == 0.0) CHECK_NOTHROW( mixture.neutral().index(x.name()) );
            if(x.charge() != 0.0) CHECK_NOTHROW( mixture.charged().index(x.name()) );
            if(x.charge()  > 0.0) CHECK_NOTHROW( mixture.cations().index(x.name()) );
            if(x.charge()  < 0.0) CHECK_NOTHROW( mixture.anions().index(x.name()) );

            // Test AqueousMixture::indicesXYZ methods
            if(x.charge() == 0.0) CHECK( contains(mixture.indicesNeutral(), species.index(x.name())) );
            if(x.charge() != 0.0) CHECK( contains(mixture.indicesCharged(), species.index(x.name())) );
            if(x.charge()  > 0.0) CHECK( contains(mixture.indicesCations(), species.index(x.name())) );
            if(x.charge()  < 0.0) CHECK( contains(mixture.indicesAnions(), species.index(x.name())) );
        }

        // Test AqueousMixture::water method
        CHECK( mixture.water().formula().equivalent("H2O") );

        // Test AqueousMixture::indexWater method
        CHECK( mixture.indexWater() == species.index("H2O") );

        // Test AqueousMixture::charges method
        CHECK( mixture.charges().size() == species.size() );
        for(auto i = 0; i < species.size(); ++i)
            CHECK( mixture.charges()[i] == species[i].charge() );

        // Test AqueousMixture::dissociationMatrix method
        CHECK( mixture.dissociationMatrix().rows() == mixture.neutral().size() );
        CHECK( mixture.dissociationMatrix().cols() == mixture.charged().size() );

        // Assemble the expected dissociation matrix for the constructed aqueous mixture
        const MatrixXd M
        { //  H+     OH-    Na+    Cl-    Ca++   Mg++   HCO3-  CO3--  K+
            {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // CO2
            {  1.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // HCl
            {  0.0,   0.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // NaCl
            {  0.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // NaOH
            {  0.0,   0.0,   0.0,   2.0,   1.0,   0.0,   0.0,   0.0,   0.0 }, // CaCl2
            {  0.0,   0.0,   0.0,   2.0,   0.0,   1.0,   0.0,   0.0,   0.0 }, // MgCl2
            {  0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   1.0,   0.0 }, // CaCO3
            {  0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   1.0,   0.0 }, // MgCO3
        };

        INFO( "dissociation matrix is\n" << mixture.dissociationMatrix() << "\nbut expected is\n" << M );
        CHECK( mixture.dissociationMatrix().isApprox(M) );

        // The temperature (in K), pressure (in Pa) and mole fractions of the aqueous species
        const auto T = 345.67;
        const auto P = 123.4e+5;
        const auto x = moleFractions(species.size());

        const auto state = mixture.state(T, P, x);

        const ArrayXd z  = mixture.charges();                            // the electric charges of the species
        const ArrayXd zc = mixture.charges()(mixture.indicesCharged());  // the electric charges of charged species

        const auto Mw = mixture.water().molarMass();                  // the molar mass of water (in kg/mol)
        const auto m  = x/(Mw * x[mixture.indexWater()]);             // the molalities of the aqueous species
        const auto mc = m(mixture.indicesCharged()).matrix();         // the molalities of the charged aqueous solutes
        const auto mn = m(mixture.indicesNeutral()).matrix();         // the molalities of the neutral aqueous solutes
        const auto ms = (mc + M.transpose() * mn).array();            // the stoichiometric molalities (as array)
        const auto Ie = 0.5 * (m * z * z).sum();                      // the effective ionic strength of the solution
        const auto Is = 0.5 * (ms * zc * zc).sum();                   // the stoichiometric ionic strength of the solution

        CHECK( state.T       == T                      );
        CHECK( state.P       == P                      );
        CHECK( state.Ie      == Approx(Ie)             );
        CHECK( state.Is      == Approx(Is)             );
        CHECK( state.rho     == Approx(997.0470390177028) );
        CHECK( state.epsilon == Approx(78.2451448082024)  );

        CHECK( state.m.isApprox(m)   );
        CHECK( state.ms.isApprox(ms) );

        WHEN("When default density and dielectric constant functions are changed")
        {
            SpeciesList species("H2O H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO3-- K+ CO2 HCl NaCl NaOH CaCl2 MgCl2 CaCO3 MgCO3");

            AqueousMixture::setDefaultWaterDensityFn([](real T, real P)
            {
                auto const wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter::Liquid);
                return wtp.D;
            });

            AqueousMixture::setDefaultWaterDielectricConstantFn([](real T, real P)
            {
                auto const wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter::Liquid);
                auto const wep = waterElectroPropsJohnsonNorton(T, P, wtp);
                return wep.epsilon;
            });

            AqueousMixture mixture(species);

            auto state = mixture.state(T, P, x);

            CHECK( state.rho     == Approx(981.650989015948) );
            CHECK( state.epsilon == Approx(63.37949243846789) );

            // Reset the default density and dielectric constant functions
            AqueousMixture::resetDefaultWaterDensityFn();
            AqueousMixture::resetDefaultWaterDielectricConstantFn();

            // Create a new AqueousMixture object
            mixture = AqueousMixture(species);

            state = mixture.state(T, P, x);

            CHECK( state.rho     == Approx(997.0470390177028) );
            CHECK( state.epsilon == Approx(78.2451448082024)  );
        }
    }
}
