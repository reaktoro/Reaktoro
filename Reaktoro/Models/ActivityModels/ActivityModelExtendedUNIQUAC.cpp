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

#include "ActivityModelExtendedUNIQUAC.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Enumerate.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>

namespace Reaktoro {

auto createActivityModelExtendedUNIQUAC(SpeciesList const& species, ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModel
{
    // Create the aqueous solution and its state object
    AqueousMixture solution(species);
    AqueousMixtureState aqstate;

    // The electrical charges of the species
    const ArrayXd z = solution.charges();

    // The electrical charges of the species squared
    const ArrayXd z2 = z * z;

    // The index of water in the solution
    const Index iw = solution.indexWater();

    // The molar mass of water (in kg/mol)
    auto const Mw = species[iw].molarMass();

    // The reference temperature used for the computation of the Debye-Huckel variable A
    auto const T0 = 273.15;

    // The reference temperature used for the temperature-correction of the energetic UNIQUAC parameters uᵢⱼ
    auto const Tr = 298.15;

    // The term b ≡ aB in the Debye-Huckel equation (assumed constant in the e-UNIQUAC model)
    auto const b = 1.5;

    // The co-ordination number z of the UNIQUAC model (Thomasen, 1997, equation 4.2)
    auto const Z = 10.0;

    // The indices of the species containing Debye-Huckel contribution (the indices of the charged species)
    Indices const iDH = solution.indicesCharged();

    // The indices of the species containing the surface area parameter rᵢ
    Indices iR;

    // The indices of the species containing the volume parameter qᵢ
    Indices iQ;

    // The indices of the species containing energetic binary interaction parameters uᵢⱼ
    Indices iU;

    ArrayXr r = ArrayXr::Zero(species.size()); // TODO: Initialize properly with values from parameters
    ArrayXr q = ArrayXr::Zero(species.size()); // TODO: Initialize properly with values from parameters

    MatrixXr u0 = MatrixXr::Constant(species.size(), species.size(), 1.0e+16); // The matrix u with the energy interaction parameters in the e-UNIQUAC model
    MatrixXr uT = MatrixXr::Constant(species.size(), species.size(), 1.0e+16); // The matrix u with the energy interaction parameters in the e-UNIQUAC model

    u0.diagonal().fill(0.0);
    uT.diagonal().fill(0.0);

    ArrayXr ln_gDH(species.size());
    ArrayXr ln_gC(species.size());
    ArrayXr ln_gR(species.size());

    ArrayXr ln_gCinf(species.size());
    ArrayXr ln_gRinf(species.size());

    ArrayXr thetapsi;  // The matrix-vector product tr(ψ)θ

    for(auto const& [formula, param] : params.r)
    {
        errorif(param <= 0.0, "The surface area parameter rᵢ in the Extended UNIQUAC model cannot be zero or negative: r[", formula, "] = ", param);
        auto const ispecies = species.findWithName(formula);
        if(ispecies >= species.size())
            continue;
        iR.push_back(ispecies);
        r[ispecies] = param;
    }

    for(auto const& [formula, param] : params.q)
    {
        errorif(param <= 0.0, "The volume parameter qᵢ in the Extended UNIQUAC model cannot be zero or negative: q[", formula, "] = ", param);
        auto const ispecies = species.findWithName(formula);
        if(ispecies >= species.size())
            continue;
        iQ.push_back(ispecies);
        q[ispecies] = param;
    }

    for(auto const& [formula1, formula2, params] : params.u)
    {
        auto const ispecies1 = species.findWithName(formula1);
        auto const ispecies2 = species.findWithName(formula2);
        if(ispecies1 >= species.size() || ispecies2 >= species.size())
            continue;
        iU.push_back(ispecies1);
        iU.push_back(ispecies2);
        u0(ispecies1, ispecies2) = u0(ispecies2, ispecies1) = params[0];
        uT(ispecies1, ispecies2) = uT(ispecies2, ispecies1) = params[1];
    }

    iR = unique(iR);
    iQ = unique(iQ);
    iU = unique(iU);

    auto const irw = index(iR, iw);
    auto const iqw = index(iQ, iw);

    errorifnot(irw < iR.size(), "The surface area parameter rᵢ in the Extended UNIQUAC model is expected for H2O.");
    errorifnot(iqw < iQ.size(), "The volume parameter qᵢ in the Extended UNIQUAC model is expected for H2O.");

    r = ArrayXr(r(iR));
    q = ArrayXr(q(iQ));

    u0 = MatrixXr(u0(iQ, iQ)); // Only binary interaction parameters involving species with positive qᵢ values
    uT = MatrixXr(uT(iQ, iQ)); // Only binary interaction parameters involving species with positive qᵢ values

    ArrayXr phi; // The ϕ array in the e-UNIQUAC model
    ArrayXr theta; // The θ array in the e-UNIQUAC model

    MatrixXr u = MatrixXr::Zero(iQ.size(), iQ.size()); // The matrix u with the energy interaction parameters in the e-UNIQUAC model
    MatrixXr psi = MatrixXr::Zero(iQ.size(), iQ.size()); // The matrix ψ in the e-UNIQUAC model

    ArrayXr xr;
    ArrayXr xq;

    // Shared pointers used in `props.extra` to avoid heap memory allocation for big objects
    auto aqstateptr = std::make_shared<AqueousMixtureState>();
    auto aqsolutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        auto const& T = args.T;
        auto const& P = args.P;
        auto const& x = args.x;

        auto const RT = universalGasConstant*T;

        // Evaluate the state of the aqueous solution
        auto const& aqstate = *aqstateptr = solution.state(T, P, x);

        // The ionic strength of the solution and its square root
        auto const& I = aqstate.Ie;
        auto const& sqrtI = sqrt(I);

        // The molalities of the species and their natural log
        auto const& m = aqstate.m;
        auto const& ln_m = log(m);

        // Set the state of matter of the phase (always liquid for the e-UNIQUAC model)
        props.som = StateOfMatter::Liquid;

        // Export the aqueous solution and its state via the `extra` data member
        props.extra["AqueousMixtureState"] = aqstateptr;
        props.extra["AqueousMixture"] = aqsolutionptr;

        // The mole fraction of water and its natural log
        auto const xw = x[iw];
        auto const ln_xw = log(xw);

        // Calculate the Debye-Huckel A coefficient according to equation (6) of Thomsen (2005)
        auto const A = 1.131 + 1.335e-3*(T - T0) + 1.164e-5*(T - T0)*(T - T0);
        auto const Lambda = 1.0 + b*sqrtI;
        auto const alpha = A*sqrtI/Lambda;

        xr = x(iR);
        xq = x(iQ);

        phi = (xr * r) / sum(xr * r);
        theta = (xq * q) / sum(xq * q);

        for(auto j = 0; j < iQ.size(); ++j)
            for(auto i = j; i < iQ.size(); ++i)
                u(j, i) = u(i, j) = u0(j, i) + uT(j, i)*(T - Tr);

        for(auto j = 0; j < iQ.size(); ++j)
            for(auto i = 0; i < iQ.size(); ++i)
                psi(j, i) = exp(-(u(j, i) - u(i, i))/T);

        thetapsi = tr(psi) * theta.matrix();

        // Reset the values of all activity coefficient contributions
        ln_gDH = 0.0;
        ln_gC = 0.0;
        ln_gR = 0.0;
        ln_gCinf = 0.0;
        ln_gRinf = 0.0;

        // Calculate the Debye-Huckel activity coefficients for charged species -- see equation (8) of Thomsen (2005) or equation (4-12) of Thomsen (1997)
        ln_gDH(iDH) = -z2(iDH)*alpha;

        // Calculate the Debye-Huckel activity coefficients for water species -- see equation (9) of Thomsen (2005) or equation (4-11) of Thomsen (1997)
        ln_gDH[iw] = 2*Mw*A*(Lambda - 1/Lambda - 2*log(Lambda))/(b*b*b);

        // The Debye-Huckel contribution to the excess Gibbs energy -- see equation (4-7) of Thomsen (1997)
        auto const GRTxDH = -4*xw*Mw*A*(log(Lambda) - Lambda + 1 + 0.5*b*b*I)/(b*b*b);

        // The UNIQUAC contribution to the excess Gibbs energy -- see equations (4-2) and (4-4) of Thomsen (1997)
        real GRTxUNIQUAC = 0.0;

        // Calculate the UNIQUAC combinatorial activity coefficients for the species -- see equation (10) of Thomsen (2005) and equation (13) of Hingerl et al. (2014)
        for(auto const& [i, ispecies] : enumerate(iR)) {
            ln_gC[ispecies] = log(phi[i]/xr[i]) + 1 - phi[i]/xr[i] - 0.5*Z*q[i]*(log(phi[i]/theta[i]) + 1 - phi[i]/theta[i]);
            ln_gCinf[ispecies] = log(r[i]/r[irw]) + 1 - r[i]/r[irw] - 0.5*Z*q[i]*(log((r[i]*q[iqw])/(r[irw]*q[i])) + 1 - (r[i]*q[iqw])/(r[irw]*q[i]));
            GRTxUNIQUAC += xr[i] * (log(phi[i]/xr[i]) - 0.5*Z*q[i]*log(phi[i]/theta[i]));
        }

        // Calculate the UNIQUAC residual activity coefficients for the species -- see equation (16) of Thomsen (2005) and equation (13) of Hingerl et al. (2014)
        for(auto const& [i, ispecies] : enumerate(iQ)) {
            ln_gR[ispecies] = q[i]*(1 - log(thetapsi[i]) - (theta * ArrayXr(psi.row(i))/thetapsi).sum());
            ln_gRinf[ispecies] = q[i]*(1 - log(psi(iqw, i)) - psi(i, iqw)/psi(iqw, iqw)); // Note: Thomsen (2005) always assume ψ(w,w) = 1, but here we don't necessarily; that's why psi(iqw, iqw) is used here.
            GRTxUNIQUAC -= xq[i] * q[i] * log(thetapsi[i]);
        }

        // Set the activity coefficients of the species in the unsymmetrical convention and molal scale -- see equation (4-14) of Thomsen (1997) -- note extra ln(xw) here to convert to molal scale
        props.ln_g = ln_gDH + (ln_gC - ln_gCinf) + (ln_gR - ln_gRinf) + ln_xw;

        // Set the activity coefficient of water species in the symmetrical convention and mole-fraction scale -- see equation (4-13) of Thomsen (1997)
        props.ln_g[iw] = ln_gDH[iw] + ln_gC[iw] + ln_gR[iw];

        // Set the activities of the aqueous species (in molality scale)
        props.ln_a = props.ln_g + ln_m;

        // Set the activitiy of water species (in mole fraction scale)
        props.ln_a[iw] = props.ln_g[iw] + ln_xw;

        // Set the excess Gibbs energy of the aqueous solution (in J/mol)
        props.Gx = (GRTxDH + GRTxUNIQUAC) * RT;
    };

    return fn;
}

auto ActivityModelExtendedUNIQUAC() -> ActivityModelGenerator
{
    return ActivityModelExtendedUNIQUAC(Params::embedded("ExtendedUNIQUAC.v2024.yaml"));
}

auto ActivityModelExtendedUNIQUAC(ActivityModelParamsExtendedUNIQUAC const& params) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelExtendedUNIQUAC(species, params); };
}

auto ActivityModelExtendedUNIQUAC(Params const& params) -> ActivityModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ActivityModelParams"), "Expecting ExtendedUNIQUAC activity model parameters in given Params object, but it lacks a `ActivityModelParams` section within which another section `ExtendedUNIQUAC` should exist.");
    errorif(!data.at("ActivityModelParams").exists("ExtendedUNIQUAC"), "Expecting ExtendedUNIQUAC activity model parameters in given Params object, under the section `ExtendedUNIQUAC`.");
    errorif(!data.at("ActivityModelParams").at("ExtendedUNIQUAC").isDict(), "Expecting section `ExtendedUNIQUAC` with ExtendedUNIQUAC activity model parameters to be a dictionary.");

    ActivityModelParamsExtendedUNIQUAC euniquacparams =
        data["ActivityModelParams"]["ExtendedUNIQUAC"].as<ActivityModelParamsExtendedUNIQUAC>();

    return ActivityModelExtendedUNIQUAC(euniquacparams);
}

} // namespace Reaktoro
