// This file is part of Reaktoro (https://reaktoro.org).
//
// Reaktoro is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// Reaktoro is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ReactiveTransportSolver.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Profiling.hpp>
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/SmartEquilibriumSolver.hpp>
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Transport/ChemicalField.hpp>
#include <Reaktoro/Transport/Mesh.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

struct ReactiveTransportSolver::Impl
{
    /// The chemical system common to all degrees of freedom in the chemical field
    ChemicalSystem system;

    /// The partition of the chemical system
    Partition partition;

    /// The solver for solving the transport equations
    TransportSolver transport_solver;

    /// The options for the reactive transport calculations
    ReactiveTransportOptions options;

    /// The result information of the last reactive transport time step calculation
    ReactiveTransportResult result;

    /// The equilibrium solver using conventional Gibbs energy minimization approach
    EquilibriumSolver equilibrium_solver;

    /// The equilibrium solver using a smart on-demand learning strategy
    SmartEquilibriumSolver smart_equilibrium_solver;

    /// The kinetic solver using conventional integration of the system of ODEs
    KineticSolver kinetic_solver;

    /// The kinetic solver using a smart on-demand learning strategy
    SmartKineticSolver smart_kinetic_solver;

    /// The list of chemical output objects
    std::vector<ChemicalOutput> outputs;

    /// The amounts of elements between fluid-equilibrium species on the boundary
    Vector be_bc;

    /// The amounts of a fluid-equilibrium elements on each cell of the mesh
    Matrix bef;

    /// The amounts of a solid-equilibrium elements on each cell of the mesh
    Matrix bes;

    /// The amounts of equilibrium elements on each cell of the mesh
    Matrix be;

    /// The amounts of scaled equilibrium elements on each cell of the mesh
    Matrix be_bar;

    /// The current number of steps in the solution of the reactive transport equations
    Index steps = 0;

    /// Name of the file and folder with a status output
    std::string folder;

    /// Number of kinetic elements
    Index Nk;

    /// Construct a ReactiveTransportSolver::Impl instance with given chemical system.
    explicit Impl(const ChemicalSystem& system)
    : Impl(Partition(system))
    {
    }

    /// Construct a ReactiveTransportSolver::Impl instance with given partition of the chemical system.
    explicit Impl(const Partition& partition)
    : system(partition.system()), partition(partition), equilibrium_solver(partition), smart_equilibrium_solver(partition),
      kinetic_solver(ReactionSystem(), partition), smart_kinetic_solver(ReactionSystem(), partition)
    {
        setBoundaryState(ChemicalState(system));
    }

    Impl(const ReactionSystem& reactions, const Partition& partition)
    : system(partition.system()), partition(partition), equilibrium_solver(partition), smart_equilibrium_solver(partition),
      kinetic_solver(reactions, partition), smart_kinetic_solver(reactions, partition)
    {
        setBoundaryState(ChemicalState(system));
    }

    /// Set the options for the reactive transport calculations.
    auto setOptions(const ReactiveTransportOptions& _options) -> void
    {
        this->options = _options;

        // Set options of equilibrium solvers
        equilibrium_solver.setOptions(options.equilibrium);
        smart_equilibrium_solver.setOptions(options.smart_equilibrium);

        // Set options of kinetic solvers
        kinetic_solver.setOptions(options.kinetics);
        this->options.kinetics.use_smart_equilibrium_solver = options.use_smart_equilibrium_solver; // TODO: make sure the flag is set

        // Set options of smart kinetic solvers
        smart_kinetic_solver.setOptions(options.smart_kinetics);
        this->options.smart_kinetics.use_smart_equilibrium_solver = options.use_smart_equilibrium_solver; // TODO: make sure the flag is set
    }

    /// Initialize the mesh discretizing the computational domain for reactive transport.
    auto setMesh(const Mesh& mesh) -> void
    {
        transport_solver.setMesh(mesh);
    }

    /// Initialize the velocity of the reactive transport model.
    auto setVelocity(double val) -> void
    {
        transport_solver.setVelocity(val);
    }

    /// Initialize the diffusion of the reactive transport model.
    auto setDiffusionCoeff(double val) -> void
    {
        transport_solver.setDiffusionCoeff(val);
    }

    /// Initialize boundary conditions of the reactive transport model.
    auto setBoundaryState(const ChemicalState& state) -> void
    {
        // The indices of the fluid species in equilibrium
        const auto& ifs = partition.indicesEquilibriumFluidSpecies();

        // The indices of the equilibrium elements
        const auto& iee = partition.indicesEquilibriumElements();

        // Get the amounts of equilibrium elements considering only fluid-equilibrium species
        be_bc = state.elementAmountsInSpecies(ifs)(iee);
    }

    /// Initialize time step of the reactive transport sequential algorithm.
    auto setTimeStep(double val) -> void
    {
        transport_solver.setTimeStep(val);
    }

    /// Add the output to the reactive transport modelling.
    auto output() -> ChemicalOutput
    {
        outputs.emplace_back(ChemicalOutput(system));
        return outputs.back();
    }

    /// Initialize the reactive transport solver.
    auto initialize() -> void
    {
        // Auxiliary variables
        const Mesh& mesh = transport_solver.mesh();
        const Index num_cells = mesh.numCells();
        const Index Ee = partition.numEquilibriumElements();

        // Initialize the number of the kinetic species
        Nk = partition.numKineticSpecies();

        // Initialize amount of elements in fluid and solid phases
        bef.resize(num_cells, Ee);
        bes.resize(num_cells, Ee);
        be.resize(num_cells, Ee);

        // Initialize equilibrium solver based on the parameter
        transport_solver.setOptions(options.transport);
        transport_solver.initialize();
    }

    /// Perform one time step of a reactive transport calculation.
    auto step(ChemicalField& field) -> void
    {
        // Depending on the presence of the kinetic species, run either kinetic or pure equilibrium simulations
        if(Nk)
            stepKinetics(field);
        else
            stepEquilibrium(field);
    }

    /// Perform one time step of a reactive transport calculation with pure equilibrium calculations.
    auto stepEquilibrium(ChemicalField& field) -> ReactiveTransportResult
    {
        // Initialize the result of the reactive transport step
        result = {};

        // Auxiliary variables
        const auto& mesh = transport_solver.mesh();
        const auto& num_cells = mesh.numCells();
        const auto& Ee = partition.numEquilibriumElements();
        const auto& iee = partition.indicesEquilibriumElements();
        const auto& ifs = partition.indicesEquilibriumFluidSpecies();
        const auto& iss = partition.indicesEquilibriumSolidSpecies();

        // Current representation of chemical states in a field
        auto& states = field.states();
        // Current representation of chemical states' properties
        auto& properties = field.properties();

        // Open the the file for outputting chemical states
        for(auto output : outputs)
        {
            output.suffix("-" + std::to_string(steps));
            output.open();
        }

        //---------------------------------------------------------------------------
        // Step 1: Perform a time step transport calculation for each fluid element
        //---------------------------------------------------------------------------
        tic(TRANSPORT_STEP)

        // Collect the amounts of elements in the solid and fluid species
        for(Index icell = 0; icell < num_cells; ++icell)
        {
            bef.row(icell) = states[icell].elementAmountsInSpecies(ifs)(iee);
            bes.row(icell) = states[icell].elementAmountsInSpecies(iss)(iee);
        }

        // Left boundary condition cell
        Index icell_bc = 0;

        // Get porosity of the left boundary cell
        const auto phi_bc = states[icell_bc].properties().fluidVolume().val / states[icell_bc].properties().volume().val;

        // Ensure the result of each fluid element transport calculation can be saved
        result.transport_of_element.resize(Ee);

        // Transport the elements in the fluid species in equilibrium
        for(Index ielement = 0; ielement < Ee; ++ielement)
        {
            // Scale BC with a porosity of the boundary cell
            transport_solver.setBoundaryValue(phi_bc * be_bc[ielement]);
            transport_solver.step(bef.col(ielement));

            // Save the result of this element transport calculation.
            result.transport_of_element[ielement] = transport_solver.result();
        }

        // Sum the amounts of elements distributed among fluid and solid species
        be.noalias() = bef + bes;

        result.timing.transport = toc(TRANSPORT_STEP);

        //---------------------------------------------------------------------------
        // Step 2: Perform a time step equilibrium calculation for each cell
        //---------------------------------------------------------------------------
        tic(EQUILIBRIUM_STEP)

        if(options.use_smart_equilibrium_solver)
        {
            // Ensure the result of each cell's smart equilibrium calculation can be saved
            result.smart_equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Scaling vector of element amounts
                double be_total = sum(be.row(icell));
                be_bar.noalias() = be.row(icell) / be_total;

                // Scale the initial vector of species amounts
                states[icell].scaleSpeciesAmounts(1 / be_total);

                // Solve with a smart equilibrium solver
                smart_equilibrium_solver.solve(states[icell], T, P, be_bar);

                // Scale back the vector of species amounts
                states[icell].scaleSpeciesAmounts(be_total);

                // Save the result of this cell's smart equilibrium calculation
                result.smart_equilibrium_at_cell[icell] = smart_equilibrium_solver.result();

            }
        }
        else
        {
            // Ensure the result of each cell's equilibrium calculation can be saved.
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                const auto T = states[icell].temperature();
                const auto P = states[icell].pressure();

                // Scale vector of element amounts
                double be_total = sum(be.row(icell));
                be_bar.noalias() = be.row(icell) / be_total;

                // Scale the initial vector of species amounts
                states[icell].scaleSpeciesAmounts(1 / be_total);

                // Solve with a conventional equilibrium solver
                equilibrium_solver.solve(states[icell], T, P, be_bar);

                // Scale back the vector of species amounts
                states[icell].scaleSpeciesAmounts(be_total);

                // Save the result of this cell's smart equilibrium calculation.
                result.equilibrium_at_cell[icell] = equilibrium_solver.result();
            }
        }

        result.timing.equilibrium = toc(EQUILIBRIUM_STEP);

        // Update the output files with the chemical state of every cell
        for(Index icell = 0; icell < num_cells; ++icell)
            for(auto output : outputs)
            {
                output.update(states[icell], icell);
                //output.update(states[icell], properties[icell], icell); // Less time-costly update without properties evaluation
            }

        // Output chemical states in the output files
        for(auto output : outputs)
            output.close();

        // Increment the current number of reactive transport steps
        ++steps;

        return result;
    }

    /// Perform one time step of a reactive transport calculation with kinetic/equilibrium calculations.
    auto stepKinetics(ChemicalField& field) -> ReactiveTransportResult
    {
        // Initialize the result of the reactive transport step
        result = {};

        // Auxiliary variables
        const auto& mesh = transport_solver.mesh();
        const auto& num_elements = system.numElements();
        const auto& num_cells = mesh.numCells();
        const auto& ifs = system.indicesFluidSpecies();
        const auto& iss = system.indicesSolidSpecies();
        const auto& dt = transport_solver.timeStep();
        const auto& t_start = steps * dt;

        auto& states = field.states();
        auto& properties = field.properties();

        // Open the the file for outputting chemical states
        for(auto output : outputs)
        {
            output.suffix("-" + std::to_string(steps));
            output.open();
        }

        //---------------------------------------------------------------------------
        // Step 1: Perform a time step transport calculation for each fluid element
        //---------------------------------------------------------------------------
        tic(TRANSPORT_STEP)

        // Collect the amounts of elements in the solid and fluid species
        for(Index icell = 0; icell < num_cells; ++icell)
        {
            bef.row(icell) = states[icell].elementAmountsInSpecies(ifs);
            bes.row(icell) = states[icell].elementAmountsInSpecies(iss);
        }

        // Left boundary condition cell
        Index icell_bc = 0;
        // Get porosity of the left boundary cell
        const auto phi_bc = states[icell_bc].properties().fluidVolume().val / states[icell_bc].properties().volume().val;

        // Ensure the result of each fluid element transport calculation can be saved
        result.transport_of_element.resize(num_elements);

        // Transport the elements in the fluid species
        for(Index ielement = 0; ielement < num_elements; ++ielement)
        {
            // Scale BC with a porosity of the boundary cell
            transport_solver.setBoundaryValue(phi_bc * be_bc[ielement]);
            transport_solver.step(bef.col(ielement));

            // Save the result of this element transport calculation.
            result.transport_of_element[ielement] = transport_solver.result();
        }

        // Sum the amounts of elements distributed among fluid and solid species
        be.noalias() = bef + bes;

        result.timing.transport = toc(TRANSPORT_STEP);

        //---------------------------------------------------------------------------
        // Step 2: Perform a time step kinetic calculation for each cell
        //---------------------------------------------------------------------------
        tic(KINETIC_STEP)

        if(options.use_smart_kinetic_solver)
        {
            // Ensure the result of each cell's smart kinetic calculation can be saved
            result.smart_kinetics_at_cell.resize(num_cells);
            result.smart_equilibrium_at_cell.resize(num_cells);
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                // Solve with a smart kinetic solver
                smart_kinetic_solver.solve(states[icell], t_start, dt, be.row(icell));

                // Update chemical properties of the field
                //properties[icell] = kinetic_solver.properties();

                // Save the result of this cell's smart kinetic calculation
                result.smart_kinetics_at_cell[icell] = smart_kinetic_solver.result();

                // Save the result of this cell's smart or conventional equilibrium calculation
                if(options.use_smart_equilibrium_solver)    result.smart_equilibrium_at_cell[icell] = smart_kinetic_solver.result().smart_equilibrium;
                else                                        result.equilibrium_at_cell[icell] = smart_kinetic_solver.result().equilibrium;
            }
        }
        else
        {
            // Ensure the result of each cell's kinetic calculation can be saved.
            result.kinetics_at_cell.resize(num_cells);
            result.smart_equilibrium_at_cell.resize(num_cells);
            result.equilibrium_at_cell.resize(num_cells);

            for(Index icell = 0; icell < num_cells; ++icell)
            {
                // Solve with a conventional kinetic solver
                kinetic_solver.solve(states[icell], t_start, dt, be.row(icell));

                // Update chemical properties of the field
                //properties[icell] = kinetic_solver.properties();

                // Save the result of this cell's kinetic calculation
                result.kinetics_at_cell[icell] = kinetic_solver.result();

                // Save the result of this cell's smart or conventional equilibrium calculation
                if(options.use_smart_equilibrium_solver)    result.smart_equilibrium_at_cell[icell] = kinetic_solver.result().smart_equilibrium;
                else                                        result.equilibrium_at_cell[icell] = kinetic_solver.result().equilibrium;

            }
        }

        result.timing.kinetics = toc(KINETIC_STEP);

        // Update the output files with the chemical state of every cell
        for(Index icell = 0; icell < num_cells; ++icell)
            for(auto output : outputs){
                output.update(states[icell], icell);
                //output.update(states[icell], properties[icell], icell);
                // TODO: use properties for outputting, when smart algorithm is able to do the extrapolation of the properties
            }

        // Output chemical states in the output files
        for(auto output : outputs)
            output.close();

        // Increment the current number of reactive transport steps
        ++steps;

        return result;
    }

    // Show database characteristics created by the ODML method
    auto outputSmartSolverInfo() const -> void
    {
        // Output smart kinetics solver info
        if(options.use_smart_kinetic_solver)
           smart_kinetic_solver.outputInfo();

        // Output smart equilibrium solver info
        if(options.use_smart_equilibrium_solver)
            smart_equilibrium_solver.outputInfo();
    }
};

ReactiveTransportSolver::ReactiveTransportSolver(const ChemicalSystem& system)
: pimpl(new Impl(system))
{
}

ReactiveTransportSolver::ReactiveTransportSolver(const Partition& partition)
: pimpl(new Impl(partition))
{
}

ReactiveTransportSolver::ReactiveTransportSolver(const ReactionSystem& reactions, const Partition& partition)
: pimpl(new Impl(reactions, partition))
{
}

ReactiveTransportSolver::~ReactiveTransportSolver()
{
}

auto ReactiveTransportSolver::operator=(ReactiveTransportSolver other) -> ReactiveTransportSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ReactiveTransportSolver::setOptions(const ReactiveTransportOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto ReactiveTransportSolver::setMesh(const Mesh& mesh) -> void
{
    pimpl->setMesh(mesh);
}

auto ReactiveTransportSolver::setVelocity(double val) -> void
{
    pimpl->setVelocity(val);
}

auto ReactiveTransportSolver::setDiffusionCoeff(double val) -> void
{
    pimpl->setDiffusionCoeff(val);
}

auto ReactiveTransportSolver::setBoundaryState(const ChemicalState& state) -> void
{
    pimpl->setBoundaryState(state);
}

auto ReactiveTransportSolver::setTimeStep(double val) -> void
{
    pimpl->setTimeStep(val);
}

auto ReactiveTransportSolver::output() -> ChemicalOutput
{
    return pimpl->output();
}

auto ReactiveTransportSolver::initialize() -> void
{
    pimpl->initialize();
}

auto ReactiveTransportSolver::step(ChemicalField& field) -> void
{
    return pimpl->step(field);
}

auto ReactiveTransportSolver::result() const -> const ReactiveTransportResult&
{
    return pimpl->result;
}

auto ReactiveTransportSolver::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto ReactiveTransportSolver::timeStep() const -> double
{
    return pimpl->transport_solver.timeStep();
}

auto ReactiveTransportSolver::outputSmartSolverInfo() const -> void
{
    return pimpl->outputSmartSolverInfo();
}

} // namespace Reaktoro
