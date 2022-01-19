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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSetup.hpp>
using namespace Reaktoro;

using autodiff::jacobian;
using autodiff::wrt;
using autodiff::at;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSetup", "[EquilibriumSetup]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto idx = [&](auto speciesname) { return system.species().index(speciesname); };

    const auto Wn = system.formulaMatrix();
    const auto Ne = Wn.rows();
    const auto Nn = Wn.cols();

    /// Mock indices for variables in x = (n, q) that are currently basic
    VectorXl ibasicvars = VectorXl{{0, 1, 2, 3}};

    EquilibriumSpecs specs(system);

    EquilibriumOptions options;

    ChemicalState state(system);
    state.setTemperature(50.0, "celsius");
    state.setPressure(100.0, "bar");
    state.setSpeciesAmount("H2O(aq)", 55.0, "mol");
    state.setSpeciesAmount("Na+(aq)", 0.1, "mol");
    state.setSpeciesAmount("Cl-(aq)", 0.1, "mol");
    state.setSpeciesAmount("CO2(g)", 1.0, "mol");
    state.setSpeciesAmount("O2(g)", 0.1, "mol");
    state.setSpeciesAmount("CaCO3(s)", 0.2, "mol");

    SECTION("Checking matrices `Aex` and `Aep` of the optimization problem")
    {
        WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
        {
            specs.temperature();
            specs.pressure();

            EquilibriumSetup setup(specs);

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            CHECK( Aex == Wn );
            CHECK( Aep.size() == 0 );
        }

        WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
        {
            specs.temperature();
            specs.volume();

            EquilibriumSetup setup(specs);

            const auto Np = 1;

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            CHECK( Aex == Wn );
            CHECK( Aep == zeros(Ne, Np) );
        }

        WHEN("volume and internal energy are input variables - the entropy maximization formulation")
        {
            specs.volume();
            specs.internalEnergy();

            EquilibriumSetup setup(specs);

            const auto Np = 2;

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            CHECK( Aex == Wn );
            CHECK( Aep == zeros(Ne, Np) );
        }

        WHEN("temperature, pressure, and pH are input variables")
        {
            specs.temperature();
            specs.pressure();
            specs.pH();

            EquilibriumSetup setup(specs);

            const auto Np = 0;
            const auto Nq = 1;

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            const auto Aen = Aex.leftCols(Nn);
            const auto Aeq = Aex.rightCols(Nq);

            CHECK( Aen == Wn );
            CHECK( Aeq == Wn.col(system.species().index("H+(aq)")) );
            CHECK( Aep == zeros(Ne, Np) );
        }

        WHEN("volume, entropy, and activity[CO2(g)] are input variables")
        {
            specs.volume();
            specs.entropy();
            specs.activity("CO2(g)");

            EquilibriumSetup setup(specs);

            const auto Np = 2;
            const auto Nq = 1;

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            const auto Aen = Aex.leftCols(Nn);
            const auto Aeq = Aex.rightCols(Nq);

            CHECK( Aen == Wn );
            CHECK( Aeq == Wn.col(system.species().index("CO2(g)")) );
            CHECK( Aep == zeros(Ne, Np) );
        }

        WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables")
        {
            specs.temperature();
            specs.pressure();
            specs.volume();
            specs.internalEnergy();
            specs.pH();
            specs.pE();
            specs.openTo("CO2");
            specs.openTo("CH4");

            EquilibriumSetup setup(specs);

            const auto Np = 2;
            const auto Nq = 2;

            const auto Aex = setup.assembleMatrixAex();
            const auto Aep = setup.assembleMatrixAep();

            const auto Aen = Aex.leftCols(Nn);
            const auto Aeq = Aex.rightCols(Nq);

            CHECK( Aen == Wn );
            CHECK( Aeq.col(0) == Wn.col(system.species().index("H+(aq)")) );
            CHECK( Aeq.col(1) == Wn.col(system.species().index("e-(aq)")) );
            CHECK( Aep.col(0) == Wn.col(system.species().index("CO2(g)")) );
            CHECK( Aep.col(1) == Wn.col(system.species().index("CH4(g)")) );
        }
    }

    SECTION("Checking functions to set/get options")
    {
        EquilibriumSpecs specs(system);
        specs.temperature();
        specs.pressure();

        EquilibriumSetup setup(specs);

        EquilibriumOptions options;

        CHECK( setup.options().epsilon == options.epsilon );

        options.epsilon = 1.0e-12;
        setup.setOptions(options);

        CHECK( setup.options().epsilon == options.epsilon );
    }

    SECTION("Checking functions to evaluate lower and upper bounds of the variables")
    {
        EquilibriumOptions options;
        options.epsilon = 1.0e-10;

        EquilibriumRestrictions restrictions(system);

        WHEN("no restrictions are actually imposed and there are no q control variables")
        {
            EquilibriumSpecs specs(system);
            specs.temperature();
            specs.pressure();

            EquilibriumSetup setup(specs);
            setup.setOptions(options);

            const auto xlower = setup.assembleLowerBoundsVector(restrictions, state);
            const auto xupper = setup.assembleUpperBoundsVector(restrictions, state);

            CHECK( xlower.size() == Nn );
            CHECK( xupper.size() == Nn );

            CHECK( (xlower.array() == options.epsilon).all() );
            CHECK( (xupper.array() == inf).all() );
        }

        WHEN("restrictions are imposed and there are ")
        {
            EquilibriumSpecs specs(system);
            specs.temperature();
            specs.pressure();
            specs.pH(); // this adds one q control variable for the amount of titrant H+
            specs.pE(); // this adds one q control variable for the amount of titrant e-

            const auto Nq = 2; // the number of q control variables

            EquilibriumSetup setup(specs);
            setup.setOptions(options);

            restrictions.cannotDecrease("Na+(aq)");
            restrictions.cannotDecreaseBelow("O2(g)", 1.0, "umol");
            restrictions.cannotDecreaseBelow("Cl-(aq)", 3.0, "umol");
            restrictions.cannotDecreaseBelow("H+(aq)", 1e-50, "mol");

            restrictions.cannotReact("CaCO3(s)");
            restrictions.cannotIncrease("CO2(g)");
            restrictions.cannotIncreaseAbove("O2(g)", 2.0, "mmol");
            restrictions.cannotIncreaseAbove("OH-(aq)", 1e-60, "mol");

            const auto xlower = setup.assembleLowerBoundsVector(restrictions, state);
            const auto xupper = setup.assembleUpperBoundsVector(restrictions, state);

            CHECK( xlower.size() == Nn + Nq );
            CHECK( xupper.size() == Nn + Nq );

            CHECK( xlower[idx("CaCO3(s)")] == state.speciesAmount("CaCO3(s)") );
            CHECK( xlower[idx("Na+(aq)")] == state.speciesAmount("Na+(aq)") );
            CHECK( xlower[idx("O2(g)")] == 1.0e-6 );
            CHECK( xlower[idx("Cl-(aq)")] == 3.0e-6 );
            CHECK( xlower[idx("H+(aq)")] == options.epsilon ); // should not be 1e-50, which is less than options.epsilon = 1e-10

            CHECK( xupper[idx("CaCO3(s)")] == state.speciesAmount("CaCO3(s)") );
            CHECK( xupper[idx("CO2(g)")] == state.speciesAmount("CO2(g)") );
            CHECK( xupper[idx("O2(g)")] == 0.002 );
            CHECK( xupper[idx("OH-(aq)")] == options.epsilon ); // should not be 1e-60, which is less than options.epsilon = 1e-10

            CHECK( (xlower.array() > options.epsilon).count() == 4 ); // CaCO3(s), Na+(aq), O2(g), Cl-(aq)
            CHECK( (xupper.array() < inf).count() == 4 ); // CaCO3(s), CO2(g), O2(g), OH-(aq)

            CHECK( (xlower.tail(Nq).array() < 0.0).all() ); // lower bounds of q variables are -inf
            CHECK( (xupper.tail(Nq).array() > 0.0).all() ); // upper bounds of q variables are -inf

            CHECK( xlower.tail(Nq).allFinite() == false ); // lower bounds of q variables are -inf
            CHECK( xupper.tail(Nq).allFinite() == false ); // upper bounds of q variables are -inf
        }
    }

    SECTION("Checking functions to evaluate objective and constraints")
    {
        const auto T = 0.7; // purely numerical value - not physically meaningful!
        const auto P = 1.3; // purely numerical value - not physically meaningful!

        WHEN("temperature and pressure are input variables")
        {
            EquilibriumSpecs specs(system);
            specs.temperature();
            specs.pressure();

            const auto Np = 0; // there are no *p* control variables in the problem
            const auto Nq = 0; // there are no *q* control variables in the problem

            const auto n = ArrayXr::LinSpaced(Nn, 1.0, Nn);
            const auto q = ArrayXr{};
            const auto p = ArrayXr{};

            ArrayXr x(Nn + Nq); // the vector x = (n, q)
            x << n, q;

            EquilibriumSetup setup(specs);

            VectorXr w{{T, P}};

            ChemicalProps props(system);
            props.update(T, P, n);

            const auto RT = universalGasConstant * T;

            const auto u = props.speciesChemicalPotentials();

            const auto G = (n * u).sum();

            const auto tau = setup.options().epsilon * setup.options().logarithm_barrier_factor;

            const auto f  = G/RT - tau * n.log().sum();
            const VectorXd fn = u/RT - tau/n;

            setup.update(x, p, w);
            setup.updateGradX(ibasicvars);
            setup.updateGradP();

            CHECK( f  == Approx(setup.getGibbsEnergy()) );
            CHECK( fn.isApprox(setup.getGibbsGradX()) );

            CHECK( setup.getConstraintResidualsGradX().size() == 0 );
            CHECK( setup.getConstraintResidualsGradP().size() == 0 );
        }

        WHEN("temperature and pressure are not input variables")
        {
            EquilibriumSpecs specs(system);
            specs.volume();
            specs.internalEnergy();
            specs.enthalpy();
            specs.pH();
            specs.pE();
            specs.openTo("CO2");

            const auto Np = 3; // p control variables: T, P, and amount of titrant [CO2] from open to CO2 condition
            const auto Nq = 2; // q control variables: the amount of titrant [H+] from pH constraint, and amount of titrant [e-] from pE constraint
            const auto Nx = Nn + Nq;

            const auto tHp  = 1.3; // current amount of titrant [H+]
            const auto tem  = 1.9; // current amount of titrant [e-]
            const auto tCO2 = 2.7; // current amount of titrant [CO2]

            auto n = ArrayXr::LinSpaced(Nn, 1.0, Nn);
            auto q = ArrayXr{{tHp, tem}};
            auto p = ArrayXr{{T, P, tCO2}};

            ArrayXr x(Nn + Nq); // the vector x = (n, q)
            x << n, q;

            EquilibriumSetup setup(specs);

            const real V = 1.0;
            const real U = 2.0;
            const real H = 3.0;
            const real pH = 4.0;
            const real pE = 5.0;

            VectorXr w{{V, U, H, pH, pE}};

            ChemicalProps props(system);
            props.update(T, P, n);

            const auto RT = universalGasConstant * T;

            const auto u = props.speciesChemicalPotentials();

            const auto G = (n * u).sum();

            const auto tau = setup.options().epsilon * setup.options().logarithm_barrier_factor;

            setup.update(x, p, w);
            setup.updateGradX(ibasicvars);
            setup.updateGradP();

            //-------------------------------------------------------------------------------
            // Check the value of the objective function at current conditions of [n, p, q]
            //-------------------------------------------------------------------------------

            const real f = G/RT - tau * n.log().sum();
            CHECK( f == Approx(setup.getGibbsEnergy()) );

            //------------------------------------------------------------------------------------------------------------
            // Check the gradient of the objective function at current conditions of [n, p, q] with respect to x = (n, q)
            //------------------------------------------------------------------------------------------------------------
            VectorXd fx(Nx);
            auto fn = fx.head(Nn);
            auto fq = fx.tail(Nq);

            fn = u/RT - tau/n;

            auto qvars = specs.controlVariablesQ();

            fq[0] = qvars[0].fn(props, w)/RT; // the pH constraint
            fq[1] = qvars[1].fn(props, w)/RT; // the pE constraint

            CHECK( fx.isApprox(setup.getGibbsGradX()) );

            //----------------------------------------------------------------------------------------------------
            // Check the Jacobian of the gradient of the objective function at current conditions of [n, p, q]
            //----------------------------------------------------------------------------------------------------

            // Auxiliary function used to update a ChemicalProps object taking
            // decision on whether to update using exact properties or
            // approximated (using ideal thermo models). This exists so that we
            // can test if EquilibriumSetup is correctly making the necessary
            // requested simplifications in derivative calculations. For
            // example, if GibbsHessian::Exact is set in
            // EquilibriumOptions::hessian, then exact properties must be
            // computed always. If GibbsHessian::PartiallyExact is set, then
            // derivatives are exact only with respect to basic variables in n.
            auto updatePropsSpecial = [ibasicvars](auto& props, auto T, auto P, auto n, auto options)
            {
                auto iseeded = 0;
                for(; iseeded < n.size() && n[iseeded][1] != 1.0; ++iseeded);

                if(iseeded < n.size()) // check there is a seeded variable in n
                {
                    auto is_nbasic_seeded = false;
                    for(auto i : ibasicvars)  // determine if seeded variable in n is a basic variable
                        if(n[i][1] == 1.0) {
                            is_nbasic_seeded = true;
                            break;
                        }

                    switch(options.hessian)
                    {
                    case GibbsHessian::Exact:
                        props.update(T, P, n);
                        break;
                    case GibbsHessian::Approx:
                    case GibbsHessian::ApproxDiagonal:
                        props.updateIdeal(T, P, n);
                        break;
                    case GibbsHessian::PartiallyExact:
                    default:
                        if(is_nbasic_seeded)
                            props.update(T, P, n); // compute exact properties if basic seeded variable
                        else props.updateIdeal(T, P, n); // compute approx properties if seeded variable is non-basic
                    }
                }
                else props.update(T, P, n); // in case there is no seeded variable in n
            };

            auto gfn = [&system, &Nn, &Nq, &Nx, &tau, &specs, &w, &options, &updatePropsSpecial](ArrayXrConstRef x, ArrayXrConstRef p)
            {
                ChemicalProps auxprops(system);
                const auto T = p[0]; // in tested equilibrium specs, p[0] is temperature (this is not necessarily always true!)
                const auto P = p[1]; // in tested equilibrium specs, p[1] is pressure (this is not necessarily always true!)
                const auto n = x.head(Nn);

                updatePropsSpecial(auxprops, T, P, n, options);

                const auto u = auxprops.speciesChemicalPotentials();
                const auto RT = universalGasConstant * T;

                VectorXr g(Nx);
                auto gn = g.head(Nn);
                auto gq = g.tail(Nq);

                gn = u/RT - tau/n.array();

                auto qvars = specs.controlVariablesQ();

                gq[0] = qvars[0].fn(auxprops, w)/RT; // the pH constraint
                gq[1] = qvars[1].fn(auxprops, w)/RT; // the pE constraint

                return g;
            };

            SECTION("checking Hxx and Hxp derivatives when using all possible GibbsHessian modes")
            {
                GibbsHessian modes[] = {
                    GibbsHessian::Exact,
                    GibbsHessian::PartiallyExact,
                    GibbsHessian::Approx,
                    GibbsHessian::ApproxDiagonal
                };

                for(auto mode : modes)
                {
                    options.hessian = mode;

                    setup.setOptions(options);
                    setup.update(x, p, w);
                    setup.updateGradX(ibasicvars);
                    setup.updateGradP();

                    const auto Hxx = jacobian(gfn, wrt(x), at(x, p));
                    const auto Hxp = jacobian(gfn, wrt(p), at(x, p));

                    CHECK( Hxx.isApprox(setup.getGibbsHessianX()) );
                    CHECK( Hxp.isApprox(setup.getGibbsHessianP()) );
                }
            }

            //-------------------------------------------------------------------------------------------------
            // Check the value of the constraint functions of equation type at current conditions of [n, p, q]
            //-------------------------------------------------------------------------------------------------
            VectorXd v(Np);

            v[0] = specs.constraintsEquationType()[0].fn(props, w); // the volume constraint equation
            v[1] = specs.constraintsEquationType()[1].fn(props, w); // the internal energy constraint equation
            v[2] = specs.constraintsEquationType()[2].fn(props, w); // the enthalpy constraint equation

            CHECK( v.isApprox(setup.getConstraintResiduals()) );

            //----------------------------------------------------------------------------------------------------
            // Check the Jacobian of the constraint functions of equation type at current conditions of [n, p, q]
            //----------------------------------------------------------------------------------------------------
            auto vfn = [&system, &Nn, &Np, &specs, &w, &options, &updatePropsSpecial](ArrayXrConstRef x, ArrayXrConstRef p)
            {
                ChemicalProps auxprops(system);
                const auto T = p[0]; // in tested equilibrium specs, p[0] is temperature (this is not necessarily always true!)
                const auto P = p[1]; // in tested equilibrium specs, p[1] is pressure (this is not necessarily always true!)
                const auto n = x.head(Nn);

                updatePropsSpecial(auxprops, T, P, n, options);

                VectorXr v(Np);
                v[0] = specs.constraintsEquationType()[0].fn(auxprops, w); // the volume constraint equation
                v[1] = specs.constraintsEquationType()[1].fn(auxprops, w); // the internal energy constraint equation
                v[2] = specs.constraintsEquationType()[2].fn(auxprops, w); // the enthalpy constraint equation

                return v;
            };

            SECTION("checking Vpx and Vpp derivatives when using all possible GibbsHessian modes")
            {
                GibbsHessian modes[] = {
                    GibbsHessian::Exact,
                    GibbsHessian::PartiallyExact,
                    GibbsHessian::Approx,
                    GibbsHessian::ApproxDiagonal
                };

                for(auto mode : modes)
                {
                    options.hessian = mode;

                    setup.setOptions(options);
                    setup.update(x, p, w);
                    setup.updateGradX(ibasicvars);
                    setup.updateGradP();

                    const auto Vpx = jacobian(vfn, wrt(x), at(x, p));
                    const auto Vpp = jacobian(vfn, wrt(p), at(x, p));

                    CHECK( Vpx.isApprox(setup.getConstraintResidualsGradX()) );
                    CHECK( Vpp.isApprox(setup.getConstraintResidualsGradP()) );
                }
            }
        }
    }
}
