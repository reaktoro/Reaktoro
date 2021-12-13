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

#ifndef REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP
#define REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP

// C++ includes
#include <map>
#include <memory>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Models/PhaseChemicalModel.hpp>

namespace Reaktoro {

// Forward declarations
class AqueousMixture;
class EUNIQUACParams;

/// Return an equation of state for an aqueous phase based on the E-UNIQUAC activity model.
/// @param mixture The aqueous mixture instance
/// @param params The parameters for the E-UNIQUAC activity model.
/// @return The activity model function for the aqueous phase
/// @see AqueousMixture, EUNIQUACParams, PhaseChemicalModel
auto aqueousChemicalModelEUNIQUAC(const AqueousMixture& mixture, const EUNIQUACParams& params) -> PhaseChemicalModel;

class EUNIQUACParams
{
public:
    /// Construct a default EUNIQUACParams instance.
    EUNIQUACParams();

    /// Get UNIQUAC r_i parameter given an species name
    auto ri(const std::string& name) const -> double;

    auto ri(const std::string& name, double value) -> void;

    auto ri(const std::map<std::string, double>& pairs) -> void;

    auto ri() const -> std::map<std::string, double>;

    /// Get UNIQUAC q_i parameter given an species name
    auto qi(const std::string& name) const -> double;

    auto qi(const std::string& name, double value) -> void;

    auto qi(const std::map<std::string, double>& pairs) -> void;

    auto qi() const -> std::map<std::string, double>;

    /// Get zeroth order energetic BIP values (Uij_0)
    auto uij_0(const std::string& first_species_name, const std::string& second_species_name) const -> double;

    auto uij_0() const -> MatrixXd;

    // TODO: enable setters
//    /// Set zeroth order energetic BIP values (Uij_0)
//    auto uij_0(const std::map<std::string, std::string>& species_pair_names, double value) -> void;

    /// Get first order energetic BIP values (Uij_T)
    auto uij_T(const std::string& first_species_name, const std::string& second_species_name) const -> double;

    auto uij_T() const -> MatrixXd;

    // TODO: enable setters
//    /// Set first order energetic BIP values (Uij_T)
//    auto uij_T(const std::map<std::string, std::string>& species_pair_names, double value) -> void;

    /// Set both uij_0 and uij_T given an species id map.
    auto set_uij_bips(
        const MatrixXd& uij_0_values,
        const MatrixXd& uij_T_values,
        const std::map<std::string, int>& species_id_map) -> void;

    /// Set E-UNIQUAC parameters values according to DTU values
    auto setDTUvalues() -> void;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

}

#endif //REAKTORO_AQUEOUSCHEMICALMODELEUNIQUAC_HPP
