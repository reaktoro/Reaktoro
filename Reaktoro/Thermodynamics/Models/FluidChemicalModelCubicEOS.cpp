// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "GaseousChemicalModelCubicEOS.hpp"

// C++ includes
#include <map>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Thermodynamics/EOS/CubicEOS.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>

namespace Reaktoro {
    namespace {

        auto fluidChemicalModelCubicEOS(const FluidMixture& mixture, CubicEOS::Model modeltype) -> PhaseChemicalModel
        {
            // The number of gases in the mixture
            const unsigned nspecies = mixture.numSpecies();

            // Get the the critical temperatures, pressures and acentric factors of the gases
            std::vector<double> Tc, Pc, omega;
            for (FluidSpecies species : mixture.species())
            {
                Tc.push_back(species.criticalTemperature());
                Pc.push_back(species.criticalPressure());
                omega.push_back(species.acentricFactor());
            }

            // Initialize the CubicEOS instance
            CubicEOS eos(nspecies);
            if (mixture.species()[0].name().find("(liq)") != std::string::npos)
                eos.setPhaseAsLiquid();
            else
                eos.setPhaseAsVapor();
            eos.setCriticalTemperatures(Tc);
            eos.setCriticalPressures(Pc);
            eos.setAcentricFactors(omega);
            eos.setModel(modeltype);

            // The state of the gaseous mixture
            FluidMixtureState state;

            // Define the chemical model function of the gaseous phase
            PhaseChemicalModel model = [=](PhaseChemicalModelResult& res, Temperature T, Pressure P, VectorConstRef n) mutable
            {
                // Evaluate the state of the gaseous mixture
                state = mixture.state(T, P, n);

                // The mole fractions of the species
                const auto& x = state.x;

                // Evaluate the CubicEOS object function
                const CubicEOS::Result eosres = eos(T, P, x);

                // The ln of mole fractions
                const ChemicalVector ln_x = log(x);

                // The ln of pressure in bar units
                const ThermoScalar ln_Pbar = log(1e-5 * P);

                // Create an alias to the ln fugacity coefficients
                const auto& ln_phi = eosres.ln_fugacity_coefficients;

                // Fill the chemical properties of the gaseous phase
                res.ln_activity_coefficients = ln_phi;
                res.ln_activities = ln_phi + ln_x + ln_Pbar;
                res.molar_volume = eosres.molar_volume;
                res.residual_molar_gibbs_energy = eosres.residual_molar_gibbs_energy;
                res.residual_molar_enthalpy = eosres.residual_molar_enthalpy;
                res.residual_molar_heat_capacity_cp = eosres.residual_molar_heat_capacity_cp;
                res.residual_molar_heat_capacity_cv = eosres.residual_molar_heat_capacity_cv;
            };

            return model;
        }

    } // namespace

    auto fluidChemicalModelVanDerWaals(const FluidMixture& mixture) -> PhaseChemicalModel
    {
        return fluidChemicalModelCubicEOS(mixture, CubicEOS::VanDerWaals);
    }

    auto fluidChemicalModelRedlichKwong(const FluidMixture& mixture) -> PhaseChemicalModel
    {
        return fluidChemicalModelCubicEOS(mixture, CubicEOS::RedlichKwong);
    }

    auto fluidChemicalModelSoaveRedlichKwong(const FluidMixture& mixture) -> PhaseChemicalModel
    {
        return fluidChemicalModelCubicEOS(mixture, CubicEOS::SoaveRedlichKwong);
    }

    auto fluidChemicalModelPengRobinson(const FluidMixture& mixture) -> PhaseChemicalModel
    {
        return fluidChemicalModelCubicEOS(mixture, CubicEOS::PengRobinson);
    }

} // namespace Reaktoro
