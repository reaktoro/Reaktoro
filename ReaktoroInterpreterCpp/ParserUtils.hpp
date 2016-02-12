// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <istream>
#include <string>
#include <tuple>
#include <vector>

// yaml-cpp includes
#include <yaml-cpp/yaml.h>

namespace Reaktoro {

/// A type used to represent a pair (value, units).
struct ValueUnits
{
    /// The value
    double value;

    /// The units
    std::string units;
};

/// A type used to represent a triplet (entity, value, units).
struct EntityValueUnits
{
    /// The value
    double value;

    /// The units
    std::string units;

    /// The entity name
    std::string entity;
};

/// A type used to represent a triplet (value, units, entity) used in a Mixture node.
struct MixtureCompound : EntityValueUnits
{
    MixtureCompound();
    MixtureCompound(std::string str);
};

/// A type used to represent a mixture node, with description of mixture compounds.
using MixtureNode = std::vector<MixtureCompound>;

/// A type used to represent an equilibrium constraint.
struct EquilibriumConstraintNode
{
    struct Base
    {
        /// The value of the constraint.
        double value;

        /// The units of the value (if any)
        std::string units;

        /// The entity name to which the constraint is applied (if any)
        std::string entity;

        /// The titrant name used to control the constraint (if any).
        std::string titrant1;

        /// The additional titrant name used to control the constraint (if any).
        /// If titrant1 and titrant2 are non-empty, then they are mutually exclusive,
        /// which means only one will exist in positive amounts, while the other is zero.
        std::string titrant2;
    };

    struct pH              : Base {};
    struct SpeciesAmount   : Base {};
    struct SpeciesActivity : Base {};
    struct PhaseAmount     : Base {};
    struct PhaseVolume     : Base {};
};

/// A type used to represent a plot.
struct PlotNode
{
    /// The name of the plot file.
    std::string name;

    /// The quantity to be plot along the x-axis.
    std::string x;

    /// The quantities to be plot along the y-axis.
    std::string y;

    /// The label used for the x-axis.
    std::string xlabel;

    /// The label used for the y-axis.
    std::string ylabel;

    /// The titles of the data plot along the y-axis.
    std::string ytitles;

    /// The settings for the key.
    std::string key;
};

/// A type used to represent an equilibrium calculation.
struct EquilibriumNode
{
    /// The name of the chemical state where this equilibrium calculation is saved.
    std::string stateid = "State";

    /// The temperature for the equilibrium calculation.
    ValueUnits temperature = {25.0, "celsius"};

    /// The pressure for the equilibrium calculation.
    ValueUnits pressure = {1.0, "bar"};

    /// The mixture definition for the equilibrium calculation.
    MixtureNode mixture;

    /// The pH constraints (only the last one used)
    std::vector<EquilibriumConstraintNode::pH> pH;

    /// The species amount constraints
    std::vector<EquilibriumConstraintNode::SpeciesAmount> species_amounts;

    /// The species activity constraints
    std::vector<EquilibriumConstraintNode::SpeciesActivity> species_activities;

    /// The phase amount constraints
    std::vector<EquilibriumConstraintNode::PhaseAmount> phase_amounts;

    /// The phase volume constraints
    std::vector<EquilibriumConstraintNode::PhaseVolume> phase_volumes;

    /// The list of inert species in the equilibrium calculation.
    std::vector<EntityValueUnits> inert_species;

    /// The list of inert phases in the equilibrium calculation.
    std::vector<std::string> inert_phases;
};

/// A type used to represent a kinetic calculation.
struct KineticsNode
{
    /// The name of the chemical state where this kinetic calculation is saved.
    std::string stateid = "State";

    /// The name of the initial chemical state from where this kinetic calculation should start.
    std::string initial_condition = "State";

    /// The names of the species that are controlled by kinetics
    std::vector<std::string> kinetic_species;

    /// The duration of the kinetic calculation
    ValueUnits duration;

    /// The plots to be executed during the calculation.
    std::vector<PlotNode> plots;
};

/// A type used to represent a mineral reaction.
struct MineralReactionNode
{
    /// The name of the mineral for this reaction.
    std::string mineral;

    /// The equation of this mineral reaction.
    std::string equation;

    /// The kinetic mechanisms of this mineral reaction.
    std::vector<std::string> mechanisms;

    /// The specific surface area of this mineral.
    ValueUnits ssa;
};

/// Return a preprocessed input script that conforms with YAML rules.
auto preprocess(std::string script) -> std::string;

/// Return a preprocessed input script that conforms with YAML rules.
auto preprocess(std::istream& stream) -> std::string;

/// The alias type to YAML::Node.
using Node = YAML::Node;

/// The type used to represent a node processor function.
using ProcessFunction = std::function<void(const Node&)>;

auto str(const Node& node) -> std::string;
auto keynode(const Node& node) -> Node;
auto valnode(const Node& node) -> Node;
auto keyword(const Node& node) -> std::string;
auto identifier(const Node& node) -> std::string;

/// Return a joined string with a node string representation.
auto operator+(std::string str, const Node& node) -> std::string;

/// Return a joined string with a node string representation.
auto operator+(const Node& node, std::string str) -> std::string;

auto operator>>(const Node& node, ValueUnits& x) -> void;
auto operator>>(const Node& node, EntityValueUnits& x) -> void;
auto operator>>(const Node& node, MixtureCompound& x) -> void;
auto operator>>(const Node& node, MixtureNode& x) -> void;
auto operator>>(const Node& node, EquilibriumConstraintNode::pH& x) -> void;
auto operator>>(const Node& node, EquilibriumConstraintNode::SpeciesAmount& x) -> void;
auto operator>>(const Node& node, EquilibriumConstraintNode::SpeciesActivity& x) -> void;
auto operator>>(const Node& node, EquilibriumConstraintNode::PhaseAmount& x) -> void;
auto operator>>(const Node& node, EquilibriumConstraintNode::PhaseVolume& x) -> void;
auto operator>>(const Node& node, PlotNode& x) -> void;
auto operator>>(const Node& node, EquilibriumNode& x) -> void;
auto operator>>(const Node& node, KineticsNode& x) -> void;
auto operator>>(const Node& node, MineralReactionNode& x) -> void;

} // namespace Reaktoro
