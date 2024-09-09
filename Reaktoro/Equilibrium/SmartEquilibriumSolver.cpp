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

#include "SmartEquilibriumSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>

namespace Reaktoro {
namespace detail {

/// Compute the step-round value of a given number.
/// This method computes the step-round value of a given number for a given step length separating the set of step-rounded values.
/// For example, if `step = 5` and `7.5 <= num <= 12.5`, `sround(num, step) => 10`.
/// @param num The number for which its step-rounded value is sought
/// @param step The length between step-rounded values
auto sround(double num, double step) -> long
{
    return round(num / step) * step;
}

} // namespace detail

struct SmartEquilibriumSolver::Impl
{
    EquilibriumSolver solver;

    EquilibriumSensitivity sensitivity;

    EquilibriumConditions conditions;

    SmartEquilibriumOptions options;

    SmartEquilibriumResult result;

    /// The temperature-pressure grid containing learned calculations for speficic temperature-pressure intervals.
    SmartEquilibriumSolver::Grid grid;

    /// Construct a SmartEquilibriumSolver::Impl object with given equilibrium problem specifications.
    Impl(EquilibriumSpecs const& specs)
    : solver(specs), sensitivity(specs), conditions(specs)
    {
        // Initialize the equilibrium solver with the default options
        setOptions(options);
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS
    //
    //=================================================================================================================

    auto solve(ChemicalState& state) -> SmartEquilibriumResult
    {
        conditions.temperature(state.temperature());
        conditions.pressure(state.pressure());
        return solve(state, conditions);
    }

    auto solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
    {
        tic(SOLVE_STEP)

        // Save a backup state in case the smart prediction fails.
        const auto statebkp = state;

        // Reset the result of the last smart equilibrium calculation
        result = {};

        // Perform a smart prediction of the chemical state
        timeit( predict(state, conditions), result.timing.prediction= )

        // Perform a learning step if the smart prediction is not satisfactory
        if (!result.prediction.accepted) {
            state = statebkp;
            timeit(learn(state, conditions), result.timing.learning = )
        }

        result.timing.solve = toc(SOLVE_STEP);

        return result;
    }

    auto solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // CHEMICAL EQUILIBRIUM METHODS WITH SENSITIVITY CALCULATION
    //
    //=================================================================================================================

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    auto solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
    {
        errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
        return {};
    }

    //=================================================================================================================
    //
    // LEARN AND PREDICT METHODS
    //
    //=================================================================================================================

    /// Perform a learning operation in which a full chemical equilibrium calculation is performed.
    auto learn(ChemicalState& state, EquilibriumConditions const& conditions) -> void
    {
        //---------------------------------------------------------------------
        // GIBBS ENERGY MINIMIZATION CALCULATION DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP)

        // Perform a full chemical equilibrium solve with sensitivity derivatives calculation
        result.learning.solve = solver.solve(state, sensitivity, conditions);

        result.timing.learning_solve = toc(EQUILIBRIUM_STEP);

        // Store a predictor only if chemical equilibrium succeded
        if (!result.learning.solve.succeeded()) {
            return;
        }

        //---------------------------------------------------------------------
        // STORAGE STEP DURING THE LEARNING PROCESS
        //---------------------------------------------------------------------
        tic(STORAGE_STEP)

        // Create an equilibrium predictor object with computed equilibrium state and its sensitivities
        EquilibriumPredictor predictor(state, sensitivity);

        // Round temperature and pressure according to their respective step lengths for discretization
        const auto iT = detail::sround(state.temperature().val(), options.temperature_step);
        const auto iP = detail::sround(state.pressure().val(), options.pressure_step);

        // Get a mutable reference to an existing temperature-pressure cell or create a new one
        auto& cell = grid.cells[{iT, iP}];

        // Generate the hash number for indices of primary species in the state
        const auto iprimary = state.equilibrium().indicesPrimarySpecies();
        const auto label = hashVector(iprimary);

        // Find the index of the cluster within the temperature-pressure grid cell that has the same primary species
        auto icluster = indexfn(cell.clusters, RKT_LAMBDA(cluster, cluster.label == label));

        // If cluster is found, store the new record in it, otherwise, create a new cluster
        if (icluster < cell.clusters.size())
        {
            auto& cluster = cell.clusters[icluster];
            cluster.records.push_back({ state, conditions, sensitivity, predictor });
            cluster.priority.extend();
        }
        else
        {
            // Create a new cluster within the current temperature-pressure grid cell
            Cluster cluster;
            cluster.iprimary = iprimary;
            cluster.label = label;
            cluster.records.push_back({ state, conditions, sensitivity, predictor });
            cluster.priority.extend();

            // Append the new cluster and initialize its connectivity and priority
            cell.clusters.push_back(cluster);
            cell.connectivity.extend();
            cell.priority.extend();
        }

        result.timing.learning_storage = toc(STORAGE_STEP);
    }

    /// Perform a prediction operation in which a chemical equilibrium state is predicted using a first-order Taylor approximation.
    auto predict(ChemicalState& state, EquilibriumConditions const& conditions) -> void
    {
        // Set the prediction status to false at the beginning
        result.prediction.accepted = false;

        // Skip prediction operation if no learning data exists yet
        if(grid.cells.empty())
            return;

        // Round temperature and pressure according to their respective step lengths for discretization
        const auto iT = detail::sround(state.temperature().val(), options.temperature_step);
        const auto iP = detail::sround(state.pressure().val(), options.pressure_step);

        // Find an existing temperature-pressure grid cell within which the state temperature/pressure are located
        auto it = grid.cells.find({iT, iP});

        // Skip prediction operation if no temperature-pressure grid cell with learning data exists yet
        if(it == grid.cells.end())
            return;

        // Get a mutable reference to the found temperature-pressure cell
        auto& cell = it->second;

        const auto wvals = conditions.inputValuesGetOrCompute(state);
        const auto cvals = conditions.initialComponentAmountsGetOrCompute(state);

        const auto w = wvals.cast<double>();
        const auto c = cvals.cast<double>();

        // Auxiliary vectors used in the lambda function below to avoid repeated memory allocation
        VectorXd dw;
        VectorXd dc;

        // The function that checks if a record in the grid pass the error test.
        auto pass_error_test = [&](Record const& record) mutable -> bool
        {
            // The primary species at the reference chemical state
            auto const& iprimary0 = record.state.equilibrium().indicesPrimarySpecies();

            // The equilibrium predictor calculator at the reference state
            auto const& predictor0 = record.predictor;

            const auto w0 = record.state.equilibrium().w();
            const auto c0 = record.state.equilibrium().c();

            dw = w - w0;
            dc = c - c0;

            using std::abs;

            for(auto ispecies : iprimary0)
            {
                const auto mu0 = predictor0.speciesChemicalPotentialReference(ispecies);
                const auto mu1 = predictor0.speciesChemicalPotentialPredicted(ispecies, dw, dc);
                if(abs(mu1 - mu0) >= options.reltol * abs(mu0) + options.abstol || isnan(mu0) || isnan(mu1))
                    return false;
            }

            return true;
        };

        // Generate the hash number for indices of primary species in the state
        const auto iprimary = state.equilibrium().indicesPrimarySpecies();
        const auto label = hashVector(iprimary);

        // The function that identifies the starting cluster index
        auto index_starting_cluster = [&]() -> Index
        {
            // If no primary species, then return number of clusters to trigger use of total usage counts of clusters
            if(iprimary.size() == 0)
                return cell.clusters.size();

            // Find the index of the cluster with the same set of primary species (search those with highest count first)
            for(auto icluster : cell.priority.order())
                if(cell.clusters[icluster].label == label)
                    return icluster;

            // In no cluster with the same set of primary species if found, then return number of clusters
            return cell.clusters.size();
        };

        // The index of the starting cluster
        const auto icluster = index_starting_cluster();

        // The ordering of the clusters to look for (starting with icluster)
        auto const& clusters_ordering = cell.connectivity.order(icluster);

        //---------------------------------------------------------------------
        // SEARCH STEP DURING THE PREDICTION PROCESS
        //---------------------------------------------------------------------
        tic(SEARCH_STEP)

        // Iterate over all clusters (starting with icluster)
        for(auto jcluster : clusters_ordering)
        {
            // Fetch records from the cluster and the order they have to be processed in
            auto const& records = cell.clusters[jcluster].records;
            auto const& records_ordering = cell.clusters[jcluster].priority.order();

            // Iterate over all records in current cluster (using the order based on the priorities)
            for(auto irecord : records_ordering)
            {
                auto const& record = records[irecord];

                //---------------------------------------------------------------------
                // ERROR CONTROL STEP DURING THE PREDICTION PROCESS
                //---------------------------------------------------------------------
                tic(ERROR_CONTROL_STEP)

                // Check if the current record passes the error test
                const auto success = pass_error_test(record);

                result.timing.prediction_error_control += toc(ERROR_CONTROL_STEP);

                if(success)
                {
                    //---------------------------------------------------------------------
                    // TAYLOR PREDICTION STEP DURING THE PREDICTION PROCESS
                    //---------------------------------------------------------------------
                    tic(TAYLOR_STEP)

                    auto const& predictor0 = record.predictor;

                    predictor0.predict(state, conditions);

                    result.timing.prediction_taylor = toc(TAYLOR_STEP);

                    // Check if all projected species amounts are positive or at least very small negative values
                    auto const& n = state.speciesAmounts();

                    const double nmin = n.minCoeff();
                    const double nsum = n.sum();

                    if(nmin <= options.reltol_negative_amounts * nsum)
                        continue; // continue searching for a another record that produces positive amounts only or tolerable negative values

                    // Check if projected species amounts conserve mass of chemical elements and charge within tolerance limits
                    const auto bnew = state.componentAmounts();
                    const auto bold = c.head(bnew.size());
                    const double bsum = bold.sum();
                    const double bdiffmax = (bnew - bold).cwiseAbs().maxCoeff();

                    if(bdiffmax > options.reltol_component_amount_conservation * bsum)
                        continue; // continue searching for a another record that produces mass conservation within tolerance limits

                    result.timing.prediction_search = toc(SEARCH_STEP);

                    //---------------------------------------------------------------------
                    // After the search is finished successfully
                    //---------------------------------------------------------------------

                    // Assign small positive values to all negative amounts
                    for(auto i = 0; i < n.size(); ++i)
                        if(n[i] < 0.0)
                            state.setSpeciesAmount(i, options.learning.epsilon);

                    //---------------------------------------------------------------------
                    // DATABASE PRIORITY UPDATE STEP DURING THE PREDICTION PROCESS
                    //---------------------------------------------------------------------
                    tic(PRIORITY_UPDATE_STEP)

                    // Increment priority of the current record (irecord) in the current cluster (jcluster)
                    cell.clusters[jcluster].priority.increment(irecord);

                    // Increment priority of the current cluster (jcluster) with respect to starting cluster (icluster)
                    cell.connectivity.increment(icluster, jcluster);

                    // Increment priority of the current cluster (jcluster)
                    cell.priority.increment(jcluster);

                    // Mark the predicted state as accepted
                    result.prediction.accepted = true;

                    result.timing.prediction_priority_update = toc(PRIORITY_UPDATE_STEP);

                    return;
                }
            }
        }

        result.prediction.accepted = false;
    }

    //=================================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================================

    /// Set the options of the smart equilibrium solver
    auto setOptions(SmartEquilibriumOptions const& opts) -> void
    {
        options = opts;
        solver.setOptions(opts.learning);
    }
};

SmartEquilibriumSolver::SmartEquilibriumSolver(ChemicalSystem const& system)
: pimpl(new Impl(EquilibriumSpecs::TP(system)))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(EquilibriumSpecs const& specs)
: pimpl(new Impl(specs))
{}

SmartEquilibriumSolver::SmartEquilibriumSolver(SmartEquilibriumSolver const& other)
: pimpl(new Impl(*other.pimpl))
{}

SmartEquilibriumSolver::~SmartEquilibriumSolver()
{}

auto SmartEquilibriumSolver::operator=(SmartEquilibriumSolver other) -> SmartEquilibriumSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto SmartEquilibriumSolver::solve(ChemicalState& state) -> SmartEquilibriumResult
{
    return pimpl->solve(state);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
{
    return pimpl->solve(state, conditions);
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumRestrictions is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::solve(ChemicalState& state, EquilibriumSensitivity& sensitivity, EquilibriumConditions const& conditions, EquilibriumRestrictions const& restrictions) -> SmartEquilibriumResult
{
    errorif(true, "SmartEquilibriumSolver::solve methods with given EquilibriumSensitivity is currently not supported.");
    return {};
}

auto SmartEquilibriumSolver::setOptions(SmartEquilibriumOptions const& options) -> void
{
    pimpl->setOptions(options);
}

} // namespace Reaktoro
