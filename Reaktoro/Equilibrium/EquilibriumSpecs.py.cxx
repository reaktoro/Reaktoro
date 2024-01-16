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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumSpecs(py::module& m)
{
    py::class_<ControlVariableQ>(m, "ControlVariableQ")
        .def(py::init<>())
        .def_readwrite("name", &ControlVariableQ::name)
        .def_readwrite("substance", &ControlVariableQ::substance)
        .def_readwrite("id", &ControlVariableQ::id)
        .def_readwrite("fn", &ControlVariableQ::fn)
        ;

    py::class_<ControlVariableP>(m, "ControlVariableP")
        .def(py::init<>())
        .def_readwrite("name", &ControlVariableP::name)
        .def_readwrite("substance", &ControlVariableP::substance)
        ;

    auto createEquationConstraintFn = [](py::function fn) -> EquationConstraintFn
    {
        const auto nargs = fn.attr("__code__").attr("co_argcount").cast<int>();
        switch(nargs)
        {
        case 1: return EquationConstraintFn([=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) { return fn(props).cast<real>(); });
        case 2: return EquationConstraintFn([=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) { return fn(props, w).cast<real>(); });
        case 3: return EquationConstraintFn([=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) { return fn(props, p, w).cast<real>(); });
        }
        errorif(true, "Expecting an equilibrium constraint function with either 1, 2, or 3 arguments, but got one with ", nargs, " arguments instead.");
        return {};
    };

    py::class_<EquationConstraintFn>(m, "EquationConstraintFn")
        .def(py::init<>())
        .def(py::init(createEquationConstraintFn))
        .def("__call__", &EquationConstraintFn::operator())
        .def("initialized", &EquationConstraintFn::initialized)
        ;

    py::implicitly_convertible<py::function, EquationConstraintFn>();

    py::class_<EquationConstraint>(m, "EquationConstraint")
        .def(py::init<>())
        .def_readwrite("id", &EquationConstraint::id)
        .def_readwrite("fn", &EquationConstraint::fn)
        ;

    // TODO: Remove this after ConstraintEquation is finally removed from C++.
    struct ConstraintEquationDeprecated {};
    py::class_<ConstraintEquationDeprecated>(m, "ConstraintEquation")
        .def(py::init([]() -> ConstraintEquationDeprecated { errorif(true, "ConstraintEquation has been renamed to EquationConstraint. Please make this change in your code."); return {}; }))
        ;

    py::class_<EquationConstraints>(m, "EquationConstraints")
        .def(py::init<>())
        .def_readwrite("ids", &EquationConstraints::ids)
        .def_readwrite("fn", &EquationConstraints::fn)
        ;

    py::class_<ReactivityConstraint>(m, "ReactivityConstraint")
        .def(py::init<>())
        .def_readwrite("id", &ReactivityConstraint::id)
        .def_readwrite("Kn", &ReactivityConstraint::Kn)
        .def_readwrite("Kp", &ReactivityConstraint::Kp)
        ;

    py::class_<ReactivityConstraints>(m, "ReactivityConstraints")
        .def(py::init<>())
        .def_readwrite("ids", &ReactivityConstraints::ids)
        .def_readwrite("Kn", &ReactivityConstraints::Kn)
        .def_readwrite("Kp", &ReactivityConstraints::Kp)
        ;

    py::class_<EquilibriumSpecs>(m, "EquilibriumSpecs")
        .def(py::init<ChemicalSystem const&>())
        .def_static("TP", &EquilibriumSpecs::TP, "Return specifications for a chemical equilbrium problem with given temperature (T) and pressure (P).")
        .def_static("HP", &EquilibriumSpecs::HP, "Return specifications for a chemical equilbrium problem with given enthalpy (H) and pressure (P).")
        .def_static("TV", &EquilibriumSpecs::TV, "Return specifications for a chemical equilbrium problem with given temperature (T) and volume (V).")
        .def_static("UV", &EquilibriumSpecs::UV, "Return specifications for a chemical equilbrium problem with given internal (U) energy and volume (V).")
        .def_static("SP", &EquilibriumSpecs::SP, "Return specifications for a chemical equilbrium problem with given entropy (S) and pressure (P).")
        .def_static("SV", &EquilibriumSpecs::SV, "Return specifications for a chemical equilbrium problem with given entropy (S) and volume (V).")
        .def("temperature", &EquilibriumSpecs::temperature, "Specify that the temperature of the system is given at chemical equilibrium.")
        .def("pressure", &EquilibriumSpecs::pressure, "Specify that the pressure of the system is given at chemical equilibrium.")
        .def("volume", &EquilibriumSpecs::volume, "Specify that the volume of the system is given at chemical equilibrium.")
        .def("internalEnergy", &EquilibriumSpecs::internalEnergy, "Specify that the internal energy of the system is given at chemical equilibrium.")
        .def("enthalpy", &EquilibriumSpecs::enthalpy, "Specify that the enthalpy of the system is given at chemical equilibrium.")
        .def("gibbsEnergy", &EquilibriumSpecs::gibbsEnergy, "Specify that the Gibbs energy of the system is given at chemical equilibrium.")
        .def("helmholtzEnergy", &EquilibriumSpecs::helmholtzEnergy, "Specify that the Helmholtz energy of the system is given at chemical equilibrium.")
        .def("entropy", &EquilibriumSpecs::entropy, "Specify that the entropy of the system is given at chemical equilibrium.")
        .def("charge", &EquilibriumSpecs::charge, "Specify that the electric charge is given at chemical equilibrium.")
        .def("elementAmount", &EquilibriumSpecs::elementAmount, "Specify that the amount of an element is given at chemical equilibrium.")
        .def("elementAmountInPhase", &EquilibriumSpecs::elementAmountInPhase, "Specify that the amount of an element in a phase is given at chemical equilibrium.")
        .def("elementMass", &EquilibriumSpecs::elementMass, "Specify that the mass of an element is given at chemical equilibrium.")
        .def("elementMassInPhase", &EquilibriumSpecs::elementMassInPhase, "Specify that the mass of an element in a phase is given at chemical equilibrium.")
        .def("phaseAmount", &EquilibriumSpecs::phaseAmount, "Specify that the amount of a phase is given at chemical equilibrium.")
        .def("phaseMass", &EquilibriumSpecs::phaseMass, "Specify that the mass of a phase is given at chemical equilibrium.")
        .def("phaseVolume", &EquilibriumSpecs::phaseVolume, "Specify that the volume of a phase is given at chemical equilibrium.")
        .def("unknownTemperature", &EquilibriumSpecs::unknownTemperature, "Specify that the temperature of the system is unknown in the chemical equilibrium calculation.")
        .def("unknownPressure", &EquilibriumSpecs::unknownPressure, "Specify that the pressure of the system is unknown in the chemical equilibrium calculation.")
        .def("chemicalPotential", &EquilibriumSpecs::chemicalPotential, "Specify that the chemical potential of a substance is given at chemical equilibrium.")
        .def("lnActivity", py::overload_cast<Species const&>(&EquilibriumSpecs::lnActivity), "Specify that the ln activity of a species is given at chemical equilibrium.")
        .def("lnActivity", py::overload_cast<String>(&EquilibriumSpecs::lnActivity), "Specify that the ln activity of a species is given at chemical equilibrium.")
        .def("lgActivity", &EquilibriumSpecs::lgActivity, "Specify that the lg activity of a species is given at chemical equilibrium.")
        .def("activity", &EquilibriumSpecs::activity, "Specify that the activity of a species is given at chemical equilibrium.")
        .def("fugacity", &EquilibriumSpecs::fugacity, "Specify that the fugacity of a gaseous species is given at chemical equilibrium.")
        .def("pH", &EquilibriumSpecs::pH, "Specify that the pH is given at chemical equilibrium.")
        .def("pMg", &EquilibriumSpecs::pMg, "Specify that pMg is given at chemical equilibrium.")
        .def("pE", &EquilibriumSpecs::pE, "Specify that pE is given at chemical equilibrium.")
        .def("Eh", &EquilibriumSpecs::Eh, "Specify that Eh is given at chemical equilibrium.")
        .def("openTo", &EquilibriumSpecs::openTo, "Specify that the chemical system is open to a substance.")
        .def("addUnknownTitrantAmount", &EquilibriumSpecs::addUnknownTitrantAmount, "Specify that the chemical system is open to a titrant substance and its amount is unknown.")
        .def("numInputs", &EquilibriumSpecs::numInputs, "Return the number of introduced input variables.")
        .def("numControlVariables", &EquilibriumSpecs::numControlVariables, "Return the number of all introduced control variables.")
        .def("numControlVariablesP", &EquilibriumSpecs::numControlVariablesP, "Return the number of introduced p control variables.")
        .def("numControlVariablesQ", &EquilibriumSpecs::numControlVariablesQ, "Return the number of introduced q control variables.")
        .def("numTitrants", &EquilibriumSpecs::numTitrants, "Return the number of all introduced explicit and implicit titrants.")
        .def("numTitrantsExplicit", &EquilibriumSpecs::numTitrantsExplicit, "Return the number of all introduced explicit titrants.")
        .def("numTitrantsImplicit", &EquilibriumSpecs::numTitrantsImplicit, "Return the number of all introduced implicit titrants.")
        .def("numEquationConstraints", &EquilibriumSpecs::numEquationConstraints, "Return the number of all introduced equation constraints.")
        .def("numReactivityConstraints", &EquilibriumSpecs::numReactivityConstraints, "Return the number of all introduced reactivity constraints.")
        .def("numConstraints", &EquilibriumSpecs::numConstraints, "Return the number of all introduced equation, reactivity, and chemical potential constraints.")
        .def("numConservativeComponents", &EquilibriumSpecs::numConstraints, "Return the number of the conservative components in the equilibrium problem specifications.")
        .def("namesInputs", &EquilibriumSpecs::namesInputs, "Return the names of the introduced input variables.")
        .def("namesControlVariables", &EquilibriumSpecs::namesControlVariables, "Return the names of all introduced control variables.")
        .def("namesControlVariablesP", &EquilibriumSpecs::namesControlVariablesP, "Return the names of introduced p control variables.")
        .def("namesControlVariablesQ", &EquilibriumSpecs::namesControlVariablesQ, "Return the names of introduced q control variables.")
        .def("namesTitrants", &EquilibriumSpecs::namesTitrants, "Return the names of all introduced explicit and implicit titrants.")
        .def("namesTitrantsExplicit", &EquilibriumSpecs::namesTitrantsExplicit, "Return the names of all introduced explicit titrants.")
        .def("namesTitrantsImplicit", &EquilibriumSpecs::namesTitrantsImplicit, "Return the names of all introduced implicit titrants.")
        .def("namesConstraints", &EquilibriumSpecs::namesConstraints, "Return the names of all introduced equation and chemical potential constraints.")
        .def("namesConservativeComponents", &EquilibriumSpecs::namesConservativeComponents, "Return the names of the conservative components in the equilibrium problem specifications.")
        .def("addControlVariableQ", &EquilibriumSpecs::addControlVariableQ, "Add a q control variable in the specification of the chemical equilibrium problem.")
        .def("addControlVariableP", &EquilibriumSpecs::addControlVariableP, "Add a p control variable in the specification of the chemical equilibrium problem.")
        .def("addConstraint", &EquilibriumSpecs::addConstraint, "Add an equation constraint to be satisfied at chemical equilibrium.")
        .def("addConstraints", &EquilibriumSpecs::addConstraints, "Add a system of equation constraints to be satisfied at chemical equilibrium.")
        .def("addReactivityConstraint", &EquilibriumSpecs::addReactivityConstraint, "Add a reactivity constraint to be satisfied at chemical equilibrium.")
        .def("addReactivityConstraints", &EquilibriumSpecs::addReactivityConstraints, "Add a system of reactivity constraints to be satisfied at chemical equilibrium.")
        .def("addInput", py::overload_cast<String const&>(&EquilibriumSpecs::addInput), "Add a new input variable for the chemical equilibrium problem with name `var`.")
        .def("system", &EquilibriumSpecs::system, "Return the chemical system associated with the equilibrium conditions.")
        .def("inputs", &EquilibriumSpecs::inputs, "Return the input variables in the chemical equilibrium specifications.")
        .def("isTemperatureUnknown", &EquilibriumSpecs::isTemperatureUnknown, "Return true if temperature is unknown in the chemical equilibrium specifications.")
        .def("isPressureUnknown", &EquilibriumSpecs::isPressureUnknown, "Return true if pressure is unknown in the chemical equilibrium specifications.")
        .def("indexTemperatureAmongInputVariables", &EquilibriumSpecs::indexTemperatureAmongInputVariables, "Return the index of temperature in the vector of w input variables if it is an input, otherwise Index(-1) if unknown.")
        .def("indexTemperatureAmongControlVariablesP", &EquilibriumSpecs::indexTemperatureAmongControlVariablesP, "Return the index of temperature in the vector of p control variables if it is unknown, otherwise Index(-1) if known.")
        .def("indexPressureAmongInputVariables", &EquilibriumSpecs::indexPressureAmongInputVariables, "Return the index of pressure in the vector of w input variables if it is an input, otherwise Index(-1) if unknown.")
        .def("indexPressureAmongControlVariablesP", &EquilibriumSpecs::indexPressureAmongControlVariablesP, "Return the index of pressure in the vector of p control variables if it is unknown, otherwise Index(-1) if known.")
        .def("indexInputVariable", &EquilibriumSpecs::indexInputVariable, "Return the index of a w input variable with given name if found, otherwise the number of w input variables.")
        .def("indexControlVariableP", &EquilibriumSpecs::indexControlVariableP, "Return the index of a p control variable with given name if found, otherwise the number of p control variables.")
        .def("indexControlVariableQ", &EquilibriumSpecs::indexControlVariableQ, "Return the index of a q control variable with given name if found, otherwise the number of q control variables.")
        .def("controlVariablesQ", &EquilibriumSpecs::controlVariablesQ, "Return the q control variables in the chemical equilibrium specifications.")
        .def("controlVariablesP", &EquilibriumSpecs::controlVariablesP, "Return the q control variables in the chemical equilibrium specifications.")
        .def("titrants", &EquilibriumSpecs::titrants, "Return the chemical formulas of the explicit and implicit titrant substances.")
        .def("titrantsExplicit", &EquilibriumSpecs::titrantsExplicit, "Return the chemical formulas of the explicit titrant substances.")
        .def("titrantsImplicit", &EquilibriumSpecs::titrantsImplicit, "Return the chemical formulas of the implicit titrant substances.")
        .def("equationConstraintsSingle", &EquilibriumSpecs::equationConstraintsSingle, "Return the specified single equation constraints to be satisfied at chemical equilibrium.")
        .def("equationConstraintsSystem", &EquilibriumSpecs::equationConstraintsSystem, "Return the specified systems of equation constraints to be satisfied at chemical equilibrium.")
        .def("reactivityConstraintsSingle", &EquilibriumSpecs::reactivityConstraintsSingle, "Return the specified single reactivity constraints to be satisfied at chemical equilibrium.")
        .def("reactivityConstraintsSystem", &EquilibriumSpecs::reactivityConstraintsSystem, "Return the specified systems of reactivity constraints to be satisfied at chemical equilibrium.")
        .def("assembleEquationConstraints", &EquilibriumSpecs::assembleEquationConstraints, "Return the complete system of equation constraints to be satisfied at chemical equilibrium.")
        .def("assembleReactivityConstraints", &EquilibriumSpecs::assembleReactivityConstraints, "Return the complete system of reactivity constraints to be satisfied at chemical equilibrium.")
        .def("assembleReactivityConstraintsMatrixKn", &EquilibriumSpecs::assembleReactivityConstraintsMatrixKn, "Assemble the coefficient matrix *Kn* of the reactivity constraints with respect to *n* species amount variables.")
        .def("assembleReactivityConstraintsMatrixKp", &EquilibriumSpecs::assembleReactivityConstraintsMatrixKp, "Assemble the coefficient matrix *Kp* of the reactivity constraints with respect to *p* control variables.")
        .def("assembleConservationMatrix", &EquilibriumSpecs::assembleConservationMatrix, "Assemble the conservation matrix for the chemical species with respect to the conservative components.")
        .def("assembleConservationMatrixN", &EquilibriumSpecs::assembleConservationMatrixN, "Assemble the conservation matrix for the *n* species amount variables with respect to the conservative components.")
        .def("assembleConservationMatrixQ", &EquilibriumSpecs::assembleConservationMatrixQ, "Assemble the conservation matrix for the *q* control variables with respect to the conservative components.")
        .def("assembleConservationMatrixP", &EquilibriumSpecs::assembleConservationMatrixP, "Assemble the conservation matrix for the *p* control variables with respect to the conservative components.")
        ;
}
