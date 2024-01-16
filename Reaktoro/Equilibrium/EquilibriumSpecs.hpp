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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/ChemicalFormula.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class ChemicalState;

/// Used to define a *q* control variable in a chemical equilibrium problem.
/// The *q* control variables are used to specify the chemical potential of a
/// species at chemical equilibrium. Its use implies that the system is open to
/// a titrant with **same** chemical formula of the substance whose chemical
/// potential needs to be constrained. A *q* control variable can also be used
/// to specify the activity of a species and thus pH.
struct ControlVariableQ
{
    /// The signature of functions that evaluate the prescribed chemical potential of a substance.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param p The control variables *p* in the chemical equilibrium calculation.
    /// @param w The input variables *w* in the chemical equilibrium calculation.
    using ChemicalPotentialFn = Fn<real(ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w)>;

    /// The unique name for this *q* control variable (required).
    String name;

    /// The chemical formula of the substance associated to this *q* control variable (required).
    ChemicalFormula substance;

    /// The unique identifier for the chemical potential constraint associated to this *q* control variable (required).
    String id;

    /// The chemical potential function associated to this *q* control variable (required).
    ChemicalPotentialFn fn;
};

/// Used to define a *p* control variable in a chemical equilibrium problem.
/// The *p* control variables are employed for designating temperature and/or
/// pressure as unknowns in the chemical equilibrium problem. They can also be
/// utilized to represent the quantities of substances, namely titrants, for
/// which the chemical system remains open.
struct ControlVariableP
{
    /// The signature of functions that evaluate the chemical potential of a species in terms of a *p* control variable.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param pk The current value of the corresponding *p* control variable during the equilibrium calculation.
    using ChemicalPotentialFn = Fn<real(ChemicalProps const& props, real const& pk)>;

    /// The unique name for this *p* control variable (required).
    String name;

    /// The chemical formula of the substance associated to this *p* control variable (optional).
    ChemicalFormula substance;
};

/// Used to define the function that evaluates the residual of an equation constraint in a chemical equilibrium problem.
struct EquationConstraintFn
{
    /// The main signature of functions that evaluate the residual of an equation constraint.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param p The control variables *p* in the chemical equilibrium calculation.
    /// @param w The input variables *w* in the chemical equilibrium calculation.
    using Func1 = Fn<real(ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w)>;

    /// An alternative signature of functions that evaluate the residual of an equation constraint.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param w The input variables *w* in the chemical equilibrium calculation.
    using Func2 = Fn<real(ChemicalProps const& props, VectorXrConstRef const& w)>;

    /// An alternative signature of functions that evaluate the residual of an equation constraint.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    using Func3 = Fn<real(ChemicalProps const& props)>;

    /// Construct a default EquationConstraintFn object.
    EquationConstraintFn() = default;

    /// Construct an EquationConstraintFn object with given equilibrium constraint function with signature @ref Func1.
    EquationConstraintFn(Func1 const& f) : _fn(f) {}

    /// Construct an EquationConstraintFn object with given equilibrium constraint function with signature @ref Func2.
    EquationConstraintFn(Func2 const& f) : _fn([=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) -> real { return f(props, w); }) {}

    /// Construct an EquationConstraintFn object with given equilibrium constraint function with signature @ref Func3.
    EquationConstraintFn(Func3 const& f) : _fn([=](ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) -> real { return f(props); }) {}

    /// Evaluate the residual of the equation constraint.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param p The control variables *p* in the chemical equilibrium calculation.
    /// @param w The input variables *w* in the chemical equilibrium calculation.
    auto operator()(ChemicalProps const& props, VectorXrConstRef const& p, VectorXrConstRef const& w) const -> real { return _fn(props, p, w); }

    /// Assign an equilibrium constraint function with signature @ref EquationConstraintFn::Func1 to this EquationConstraintFn object.
    auto operator=(Func1 const& f) { return *this = EquationConstraintFn(f); }

    /// Assign an equilibrium constraint function with signature @ref EquationConstraintFn::Func2 to this EquationConstraintFn object.
    auto operator=(Func2 const& f) { return *this = EquationConstraintFn(f); }

    /// Assign an equilibrium constraint function with signature @ref EquationConstraintFn::Func3 to this EquationConstraintFn object.
    auto operator=(Func3 const& f) { return *this = EquationConstraintFn(f); }

    /// Return `true` if this EquationConstraintFn object has been initialized with an equilibrium constraint function.
    auto initialized() const { return static_cast<bool>(_fn); }

    /// Convert this EquationConstraintFn object to a boolean value.
    operator bool() const { return initialized(); }

    /// The function that evaluates the residual of the equation constraint.
    Func1 _fn;
};

/// Used to define equation constraints in a chemical equilibrium problem.
struct EquationConstraint
{
    /// The unique identifier for this equation constraint.
    String id;

    /// The function defining the equation to be satisfied at chemical equilibrium.
    EquationConstraintFn fn;
};

/// Used to define equation constraints in a chemical equilibrium problem.
struct [[deprecated("ConstraintEquation has been renamed to EquationConstraint. Please make this change in your code.")]] ConstraintEquation : EquationConstraint
{};

/// Used to define a system of equation constraints in a chemical equilibrium problem.
struct EquationConstraints
{
    /// The signature of functions that evaluate the residual of the system of equation constraints.
    /// @param props The current chemical properties of the system in the equilibrium calculation.
    /// @param w The input variables *w* in the chemical equilibrium calculation.
    /// @param p The control variables *p* in the chemical equilibrium calculation.
    using Func = Fn<VectorXr(ChemicalProps const& props, VectorXrConstRef const& w, VectorXrConstRef const& p)>;

    /// The unique identifier for each equation constraint.
    Strings ids;

    /// The function defining the system of equations to be satisfied at chemical equilibrium.
    Func fn;
};

/// Used to define reactivity restrictions among species in the chemical
/// equilibrium calculation. This can be used, for example, to impose that a
/// reaction is inert and should not progress during the equilibration process.
struct ReactivityConstraint
{
    /// The unique identifier for this reactivity constraint.
    String id;

    /// The linear equation coefficients in the constraint corresponding to the species amounts variables *n*.
    VectorXd Kn;

    /// The linear equation coefficients in the constraint corresponding to the introduced control variables *p*.
    VectorXd Kp;
};

/// Used to define a system of reactivity restrictions among species in the chemical
/// equilibrium calculation. This can be used, for example, to impose that a collection of
/// reactions are inert and should not progress during the equilibration process.
struct ReactivityConstraints
{
    /// The unique identifiers for each reactivity constraint.
    Strings ids;

    /// The coefficient matrix of the linear reactivity constraint equations corresponding to the species amounts variables *n*.
    MatrixXd Kn;

    /// The coefficient matrix of the linear reactivity constraint equations corresponding to the introduced control variables *p*.
    MatrixXd Kp;
};

/// The class used to define conditions to be satisfied at chemical equilibrium.
///
/// ### Explicit Titrants
///
/// The *explicit titrants* are the titrants introduced with method
/// EquilibriumSpecs::openTo. In the code below, the names of the explicitly
/// introduced titrants are shown:
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.openTo("H2S"); // a titrant named [H2S] will be introduced.
/// specs.openTo("CO2"); // a titrant named [CO2] will be introduced.
/// ~~~
/// The amounts of these explicit titrants are unknown in the chemical
/// equilibrium problem and computed together with the amounts of species.
///
/// ### Implicit Titrants
///
/// The *implicit titrants* are the titrants introduced with methods:
///
/// * EquilibriumSpecs::chemicalPotential
/// * EquilibriumSpecs::lnActivity
/// * EquilibriumSpecs::lgActivity
/// * EquilibriumSpecs::activity
/// * EquilibriumSpecs::fugacity
/// * EquilibriumSpecs::pH
/// * EquilibriumSpecs::pMg
/// * EquilibriumSpecs::pE
/// * EquilibriumSpecs::Eh
///
/// In the code below, the names of the implicitly introduced titrants are shown:
///
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.lnActivity("Ca+2"); // a titrant named [Ca+2] will be introduced.
/// specs.fugacity("O2");     // a titrant named [O2] will be introduced.
/// specs.pH();               // a titrant named [H+] will be introduced.
/// specs.pMg();              // a titrant named [Mg+2] will be introduced.
/// specs.pE();               // a titrant named [e-] will be introduced.
/// specs.Eh();               // a titrant named [e-] will be introduced.
/// ~~~
///
/// The amounts of these implicit titrants are unknown in the chemical
/// equilibrium problem and computed together with the amounts of species.
///
/// ### Control Variables
///
/// The *control variables* in a chemical equilibrium problem are unknowns
/// introduced along with equilibrium constraints. These control variables can be:
///
/// * temperature,
/// * pressure, and
/// * amounts of explicit and implicit titrants.
///
/// The number of these control variables depend whether temperature and/or
/// pressure are unknown and if any titrant has been introduced, explicitly
/// or implicitly. For example, if EquilibriumSpecs::temperature is not
/// called, then temperature is unknown and computed in the chemical
/// equilibrium calculation. The same applies for pressure in case
/// EquilibriumSpecs::pressure is not called. Titrants are introduced
/// either explicitly, with method EquilibriumSpecs::openTo, or implicitly,
/// with methods:
///
/// * EquilibriumSpecs::chemicalPotential
/// * EquilibriumSpecs::lnActivity
/// * EquilibriumSpecs::lgActivity
/// * EquilibriumSpecs::activity
/// * EquilibriumSpecs::fugacity
/// * EquilibriumSpecs::pH
/// * EquilibriumSpecs::pMg
/// * EquilibriumSpecs::pE
/// * EquilibriumSpecs::Eh
///
/// The following example formulates a set of specifications for a chemical
/// equilibrium problem in which temperature and the amount of titrant
/// `[CO2]` are introduced *control variables*. Their values are not known
/// in advance; they are computed as part of the chemical equilibrium calculation.
///
/// ~~~{.cpp}
/// using namespace Reaktoro;
/// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
/// specs.temperature(); // temperature is an input variable in the chemical equilibrium problem.
/// specs.pressure();    // pressure is an input variable in the chemical equilibrium problem.
/// specs.volume();      // volume is an input variable in the chemical equilibrium problem.
/// specs.openTo("CO2"); // a titrant named [CO2] is introduced, and its amount computed at the end of the equilibrium calculation.
/// ~~~
class EquilibriumSpecs
{
public:
    /// Construct an EquilibriumSpecs object.
    explicit EquilibriumSpecs(ChemicalSystem const& system);

    //=================================================================================================
    //
    // STATIC METHODS TO CREATE PREDEFINED CHEMICAL EQUILIBRIUM SPECIFICATIONS
    //
    //=================================================================================================

    /// Return specifications for a chemical equilbrium problem with given *temperature (T)* and *pressure (P)*.
    static auto TP(ChemicalSystem const& system) -> EquilibriumSpecs;

    /// Return specifications for a chemical equilbrium problem with given *enthalpy (H)* and *pressure (P)*.
    static auto HP(ChemicalSystem const& system) -> EquilibriumSpecs;

    /// Return specifications for a chemical equilbrium problem with given *temperature (T)* and *volume (V)*.
    static auto TV(ChemicalSystem const& system) -> EquilibriumSpecs;

    /// Return specifications for a chemical equilbrium problem with given *internal (U) energy* and *volume (V)*.
    static auto UV(ChemicalSystem const& system) -> EquilibriumSpecs;

    /// Return specifications for a chemical equilbrium problem with given *entropy (S)* and *pressure (P)*.
    static auto SP(ChemicalSystem const& system) -> EquilibriumSpecs;

    /// Return specifications for a chemical equilbrium problem with given *entropy (S)* and *volume (V)*.
    static auto SV(ChemicalSystem const& system) -> EquilibriumSpecs;

    //=================================================================================================
    //
    // METHODS TO SPECIFY THERMODYNAMIC CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **temperature** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `T`. By calling this method, you are
    /// specifying that temperature is known in the equilibrium calculation (by default it is
    /// considered an unknown). Thus, it will not be considered as a control variable whose value
    /// needs to be computed as part of the equilibrium calculation.
    /// @see unknownTemperature
    auto temperature() -> void;

    /// Specify that the **pressure** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `P`. By calling this method, you are
    /// specifying that pressure is known in the equilibrium calculation (by default it is
    /// considered an unknown). Thus, it will not be considered as a control variable whose value
    /// needs to be computed as part of the equilibrium calculation.
    /// @see unknownPressure
    auto pressure() -> void;

    /// Specify that the **volume** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `V`. It also
    /// introduces an equation constraint with name `volume` to enforce a given
    /// volume value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto volume() -> void;

    /// Specify that the **internal energy** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `U`. It also
    /// introduces an equation constraint with name `internalEnergy` to enforce a given
    /// internal energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto internalEnergy() -> void;

    /// Specify that the **enthalpy** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `H`. It also
    /// introduces an equation constraint with name `enthalpy` to enforce a given
    /// enthalpy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto enthalpy() -> void;

    /// Specify that the **Gibbs energy** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `G`. It also
    /// introduces an equation constraint with name `gibbsEnergy` to enforce a given
    /// Gibbs energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto gibbsEnergy() -> void;

    /// Specify that the **Helmholtz energy** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `A`. It also
    /// introduces an equation constraint with name `helmholtzEnergy` to enforce a given
    /// Helmholtz energy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto helmholtzEnergy() -> void;

    /// Specify that the **entropy** of the system is given at chemical equilibrium.
    /// This method introduces one input variable with name `S`. It also
    /// introduces an equation constraint with name `entropy` to enforce a given
    /// entropy value for the chemical system at equilibrium.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints.
    auto entropy() -> void;

    /// Specify that the **electric charge** is given at chemical equilibrium.
    auto charge() -> void;

    /// Specify that the **amount of an element** is given at chemical equilibrium.
    /// @param element The name or index of the element in the chemical system.
    auto elementAmount(StringOrIndex const& element) -> void;

    /// Specify that the **amount of an element in a phase** is given at chemical equilibrium.
    /// @param element The name or index of the element in the chemical system.
    /// @param phase The name or index of the phase in the chemical system.
    auto elementAmountInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void;

    /// Specify that the **mass of an element** is given at chemical equilibrium.
    /// @param element The name or index of the element in the chemical system.
    auto elementMass(StringOrIndex const& element) -> void;

    /// Specify that the **mass of an element in a phase** is given at chemical equilibrium.
    /// @param element The name or index of the element in the chemical system.
    /// @param phase The name or index of the phase in the chemical system.
    auto elementMassInPhase(StringOrIndex const& element, StringOrIndex const& phase) -> void;

    /// Specify that the **amount of a phase** is given at chemical equilibrium.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseAmount(StringOrIndex const& phase) -> void;

    /// Specify that the **mass of a phase** is given at chemical equilibrium.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseMass(StringOrIndex const& phase) -> void;

    /// Specify that the **volume of a phase** is given at chemical equilibrium.
    /// @param phase The name or index of the phase in the chemical system.
    auto phaseVolume(StringOrIndex const& phase) -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY UNKNOWN INPUT CONDITIONS
    //
    //=================================================================================================

    /// Specify that the *temperature* of the system is unknown in the chemical equilibrium calculation.
    /// This method introduces one *p* control variable with name `T`. By calling this method, you
    /// are specifying that temperature is unknown and to be computed as part of the chemical
    /// equilibrium calculation.
    auto unknownTemperature() -> void;

    /// Specify that the *pressure* of the system is unknown in the chemical equilibrium calculation.
    /// This method introduces one *p* control variable with name `P`. By calling this method, you
    /// are specifying that pressure is unknown and to be computed as part of the chemical
    /// equilibrium calculation.
    auto unknownPressure() -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY CHEMICAL POTENTIAL CONSTRAINTS
    //
    //=================================================================================================

    /// Specify that the **chemical potential** of a substance is given at chemical equilibrium.
    /// This method introduces one input variable with name `u[substance]`
    /// (e.g., `u[H2O]` if @p substance is `"H2O"`). It also introduces a
    /// chemical potential constraint with same name to enforce a given
    /// chemical potential value for the substance at chemical equilibrium.
    /// This method also indicates that the chemical system is open to @p
    /// substance (e.g., the system is open to mass in/out of H@sub{2}O). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[H2O]`).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param substance The chemical formula of the substance (e.g., `H2O`, `CO2`, `H+`, `Mg+2`).
    auto chemicalPotential(String substance) -> void;

    /// Specify that the **ln activity** of a species is given at chemical equilibrium.
    /// This method introduces one input variable with name
    /// `lnActivity[speciesName]` (e.g., `lnActivity[CO2(aq)]` if @p species is
    /// a Species object with name `"CO2(aq)"`). It also introduces a chemical
    /// potential constraint with same name that is equivalent to enforcing a
    /// given value for the natural log activity of the species at chemical
    /// equilibrium. This method also indicates that the chemical system is
    /// open to the underlying substance of the species, not the species
    /// itself. For example, the system is open to mass in/out of CO@sub{2} if
    /// the Species object @p species is `CO2(aq)`, `CO2(g)` or `CO2(l)`). Thus
    /// an *implicit titrant* is introduced with name `[substance]` (e.g., `[CO2]`).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The chemical species as an Species object.
    auto lnActivity(Species const& species) -> void;

    /// Specify that the **ln activity** of a species is given at chemical equilibrium.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lnActivity(String species) -> void;

    /// Specify that the **lg activity** of a species is given at chemical equilibrium.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto lgActivity(String species) -> void;

    /// Specify that the **activity** of a species is given at chemical equilibrium.
    /// For more details, check the documentation of EquilibriumSpecs::lnActivity(Species const&).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param species The name of the chemical species as found in the database in use.
    /// @note The chemical species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a species with given name.
    auto activity(String species) -> void;

    /// Specify that the **fugacity** of a gaseous species is given at chemical equilibrium.
    /// This method introduces one input variable with name `f[gas]` (e.g.,
    /// `f[O2]` if @p gas is `"O2"`). It also introduces a chemical potential
    /// constraint with same name that is equivalent to enforcing a given value
    /// for the fugacity of the gas at chemical equilibrium. This method also
    /// indicates that the chemical system is open to @p gas (e.g., the system is
    /// open to mass in/out of O@sub{2}). Thus an *implicit titrant* is
    /// introduced with name `[substance]` (e.g., `[O2]`).
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param gas The name of the gaseous species as found in the database in use.
    /// @note The gaseous species does not need to be in the chemical system; only in the database.
    /// @warning An error will be thrown if the database does not contain a gaseous species with given name.
    auto fugacity(String gas) -> void;

    /// Specify that the pH is given at chemical equilibrium.
    /// This method introduces one input variable with name `pH`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pH at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of H@sup{+}. Thus an *implicit titrant* is introduced with name `[H+]`.
    /// The code below demonstrate the use of this method and its effect on
    /// the list of input variables, titrants and control variables.
    ///
    /// ~~~{.cpp}
    /// using namespace Reaktoro;
    /// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
    /// specs.enthalpy();                     // introduces input variable `H` and constraint `enthalpy`
    /// specs.pressure();                     // introduces input variable `P`
    /// specs.pH();                           // introduces input variable `pH`, constraint `pH`, and titrant `[H+]`
    /// print(specs.namesInputs());      // H, P, pH
    /// print(specs.namesTitrants());         // [H+]
    /// print(specs.namesControlVariables()); // T, n[H+]
    /// print(specs.namesContraints());       // enthalpy, pH
    /// ~~~
    ///
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `H+`.
    auto pH() -> void;

    /// Specify that pMg is given at chemical equilibrium.
    /// This method introduces one input variable with name `pMg`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pMg at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of Mg@sup{2+}. Thus an *implicit titrant* is introduced with name `[Mg+2]`.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @warning An error will be thrown if the system does not contain an aqueous species with formula `Mg+2`.
    auto pMg() -> void;

    /// Specify that pE is given at chemical equilibrium.
    /// This method introduces one input variable with name `pE`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for pE at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    auto pE() -> void;

    /// Specify that Eh is given at chemical equilibrium.
    /// This method introduces one input variable with name `Eh`. It also
    /// introduces a chemical potential constraint with same name that is
    /// equivalent to enforcing a given value for Eh at chemical equilibrium.
    /// This method also indicates that the chemical system is open to mass
    /// in/out of e@sup{-} (the electron substance). Thus an *implicit titrant*
    /// is introduced with name `[e-]`.
    /// @see EquilibriumSpecs::namesInputs, EquilibriumSpecs::namesConstraints, EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    auto Eh() -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY HOW THE CHEMICAL SYSTEM IS OPEN
    //
    //=================================================================================================

    /// Specify that the chemical system is open to a substance.
    ///
    /// Use this method to specify that the system is titrated with an unknown
    /// amount of a substance to be able to attain chemical equilibrium with
    /// given conditions. Its use introduces an *explicit titrant* with name
    /// `[substance]` (e.g., `H2S` is @p substance is `"H2S"`). The amount of
    /// this titrant is an unknown control variable which is computed in the
    /// chemical equilibrium calculation.
    ///
    /// The code below demonstrate the use of this method and its effect on
    /// the list of titrants and control variables.
    ///
    /// ~~~{.cpp}
    /// using namespace Reaktoro;
    /// EquilibriumSpecs specs(system); // for some ChemicalSystem object `system`
    /// specs.volume();
    /// specs.pressure();
    /// specs.openTo("H2S");
    /// specs.openTo("CO2");
    /// print(specs.namesTitrants()); // [H2S], [CO2]
    /// print(specs.namesControlVariables()); // T, [H2S], [CO2]
    /// ~~~
    ///
    /// @note The code above is for demonstration purposes only. The given
    /// specifications do not produce a valid chemical equilibrium problem!
    ///
    /// @see EquilibriumSpecs::namesTitrants, EquilibriumSpecs::namesControlVariables
    /// @param substance The chemical formula of the substance.
    auto openTo(ChemicalFormula const& substance) -> void;

    /// Specify that the chemical system is open to a given chemical state. // TODO: Implement EquilibriumSpecs::auto openTo(ChemicalState const& state).
    // auto openTo(ChemicalState const& state) -> void;

    /// Specify that the chemical system is open to a given material. // TODO: Implement EquilibriumSpecs::auto openTo(Material const& material).
    // auto openTo(Material const& material) -> void;

    //=================================================================================================
    //
    // METHODS TO SPECIFY ADDITIONAL UNKNOWNS
    //
    //=================================================================================================

    /// Specify that the chemical system is open to a titrant substance and its amount is unknown.
    auto addUnknownTitrantAmount(ChemicalFormula const& substance) -> void;

    //=================================================================================================
    //
    // METHODS TO GET THE NUMBER OF INTRODUCED CONSTRAINTS, INPUT VARIABLES, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the number of introduced input variables.
    auto numInputs() const -> Index;

    /// Return the number of all introduced control variables.
    auto numControlVariables() const -> Index;

    /// Return the number of introduced *p* control variables.
    auto numControlVariablesP() const -> Index;

    /// Return the number of introduced *q* control variables.
    auto numControlVariablesQ() const -> Index;

    /// Return the number of all introduced explicit and implicit titrants.
    auto numTitrants() const -> Index;

    /// Return the number of all introduced explicit titrants.
    auto numTitrantsExplicit() const -> Index;

    /// Return the number of all introduced implicit titrants.
    auto numTitrantsImplicit() const -> Index;

    /// Return the number of all introduced equation constraints.
    auto numEquationConstraints() const -> Index;

    /// Return the number of all introduced reactivity constraints.
    auto numReactivityConstraints() const -> Index;

    /// Return the number of all introduced equation, reactivity, and chemical potential constraints.
    auto numConstraints() const -> Index;

    /// Return the number of the conservative components in the equilibrium problem specifications.
    auto numConservativeComponents() const -> Index;

    //=================================================================================================
    //
    // METHODS TO GET THE NAMES OF INTRODUCED CONSTRAINTS, INPUT VARIABLES, AND CONTROL VARIABLES
    //
    //=================================================================================================

    /// Return the names of the introduced input variables.
    auto namesInputs() const -> Strings;

    /// Return the names of all introduced control variables.
    auto namesControlVariables() const -> Strings;

    /// Return the names of introduced *p* control variables.
    auto namesControlVariablesP() const -> Strings;

    /// Return the names of introduced *q* control variables.
    auto namesControlVariablesQ() const -> Strings;

    /// Return the names of all introduced explicit and implicit titrants.
    auto namesTitrants() const -> Strings;

    /// Return the names of all introduced explicit titrants.
    auto namesTitrantsExplicit() const -> Strings;

    /// Return the names of all introduced implicit titrants.
    auto namesTitrantsImplicit() const -> Strings;

    /// Return the names of all introduced equation and chemical potential constraints.
    auto namesConstraints() const -> Strings;

    /// Return the names of the conservative components in the equilibrium problem specifications.
    auto namesConservativeComponents() const -> Strings;

    //=================================================================================================
    //
    // METHODS TO ADD CONTROL VARIABLES, CONSTRAINTS, AND INPUT VARIABLES
    //
    //=================================================================================================

    /// Add a *q* control variable in the specification of the chemical equilibrium problem.
    auto addControlVariableQ(ControlVariableQ const& qvar) -> void;

    /// Add a *p* control variable in the specification of the chemical equilibrium problem.
    auto addControlVariableP(ControlVariableP const& pvar) -> void;

    /// Add an equation constraint to be satisfied at chemical equilibrium.
    auto addConstraint(EquationConstraint const& constraint) -> void;

    /// Add a system of equation constraints to be satisfied at chemical equilibrium.
    auto addConstraints(EquationConstraints const& constraints) -> void;

    /// Add a reactivity constraint to be satisfied at chemical equilibrium.
    auto addReactivityConstraint(ReactivityConstraint const& constraint) -> void;

    /// Add a system of reactivity constraints to be satisfied at chemical equilibrium.
    auto addReactivityConstraints(ReactivityConstraints const& constraints) -> void;

    /// Add a new input variable for the chemical equilibrium problem with name @p var.
    auto addInput(String const& var) -> Index;

    //=================================================================================================
    //
    // MISCELLANEOUS METHODS
    //
    //=================================================================================================

    /// Return the chemical system associated with the equilibrium conditions.
    auto system() const -> ChemicalSystem const&;

    /// Return the input variables in the chemical equilibrium specifications.
    auto inputs() const -> Strings const&;

    /// Return true if temperature is unknown in the chemical equilibrium specifications.
    auto isTemperatureUnknown() const -> bool;

    /// Return true if pressure is unknown in the chemical equilibrium specifications.
    auto isPressureUnknown() const -> bool;

    /// Return the index of temperature in the vector of *w* input variables if it is an input, otherwise `Index(-1)` if unknown.
    auto indexTemperatureAmongInputVariables() const -> Index;

    /// Return the index of temperature in the vector of *p* control variables if it is unknown, otherwise `Index(-1)` if known.
    auto indexTemperatureAmongControlVariablesP() const -> Index;

    /// Return the index of pressure in the vector of *w* input variables if it is an input, otherwise `Index(-1)` if unknown.
    auto indexPressureAmongInputVariables() const -> Index;

    /// Return the index of pressure in the vector of *p* control variables if it is unknown, otherwise `Index(-1)` if known.
    auto indexPressureAmongControlVariablesP() const -> Index;

    /// Return the index of a *w* input variable with given name if found, otherwise the number of *w* input variables.
    auto indexInputVariable(String const& name) const -> Index;

    /// Return the index of a *p* control variable with given name if found, otherwise the number of *p* control variables.
    auto indexControlVariableP(String const& name) const -> Index;

    /// Return the index of a *q* control variable with given name if found, otherwise the number of *q* control variables.
    auto indexControlVariableQ(String const& name) const -> Index;

    /// Return the *q* control variables in the chemical equilibrium specifications.
    auto controlVariablesQ() const -> Vec<ControlVariableQ> const&;

    /// Return the *q* control variables in the chemical equilibrium specifications.
    auto controlVariablesP() const -> Vec<ControlVariableP> const&;

    /// Return the chemical formulas of the explicit and implicit titrant substances.
    auto titrants() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the explicit titrant substances.
    auto titrantsExplicit() const -> Vec<ChemicalFormula>;

    /// Return the chemical formulas of the implicit titrant substances.
    auto titrantsImplicit() const -> Vec<ChemicalFormula>;

    /// Return the specified single equation constraints to be satisfied at chemical equilibrium.
    auto equationConstraintsSingle() const -> Vec<EquationConstraint> const&;

    /// Return the specified systems of equation constraints to be satisfied at chemical equilibrium.
    auto equationConstraintsSystem() const -> Vec<EquationConstraints> const&;

    /// Return the specified single reactivity constraints to be satisfied at chemical equilibrium.
    auto reactivityConstraintsSingle() const -> Vec<ReactivityConstraint> const&;

    /// Return the specified systems of reactivity constraints to be satisfied at chemical equilibrium.
    auto reactivityConstraintsSystem() const -> Vec<ReactivityConstraints> const&;

    //=================================================================================================
    //
    // METHODS THAT ASSEMBLE CONSTRAINTS AND MATRICES FOR CURRENT STATE OF EQUILIBRIUM SPECIFICATIONS
    //
    //=================================================================================================

    /// Assemble the complete system of equation constraints to be satisfied at chemical equilibrium.
    /// This method assembles a system of equation constraints that considers first the
    /// evaluation of the provided single equations via @ref addConstraint and then
    /// the provided systems of equations via method @ref addConstraints.
    auto assembleEquationConstraints() const -> EquationConstraints;

    /// Assemble the complete system of reactivity constraints to be satisfied at chemical equilibrium.
    /// This method assembles a system of reactivity constraints that considers first the
    /// evaluation of the provided single constraints via @ref addReactivityConstraint and then
    /// the provided systems of constraints via method @ref addReactivityConstraints.
    auto assembleReactivityConstraints() const -> ReactivityConstraints;

    /// Assemble the coefficient matrix *Kn* of the reactivity constraints with respect to *n* species amount variables.
    auto assembleReactivityConstraintsMatrixKn() const -> MatrixXd;

    /// Assemble the coefficient matrix *Kp* of the reactivity constraints with respect to *p* control variables.
    auto assembleReactivityConstraintsMatrixKp() const -> MatrixXd;

    /// Assemble the conservation matrix for the chemical species with respect to the conservative components.
    /// This method is equivalent to @ref assembleConservationMatrixN
    auto assembleConservationMatrix() const -> MatrixXd;

    /// Assemble the conservation matrix for the *n* species amount variables with respect to the conservative components.
    /// This conservative matrix of the chemical equilibrium problem is a matrix whose upper rows
    /// contains the formula matrix of the species with respect to elements and electric charge, and
    /// the lower rows contains the coefficients of the reactivity constraints (e.g., the
    /// stoichiometric matrix of the restricted reactions in the equilibrium computation). This
    /// matrix is used to compute the amounts of the conservative components in the chemical
    /// equilibrium problem, which are elements, electric charge, and the extend of the restricted
    /// reactions.
    auto assembleConservationMatrixN() const -> MatrixXd;

    /// Assemble the conservation matrix for the *q* control variables with respect to the conservative components.
    auto assembleConservationMatrixQ() const -> MatrixXd;

    /// Assemble the conservation matrix for the *p* control variables with respect to the conservative components.
    auto assembleConservationMatrixP() const -> MatrixXd;

private:
    /// The chemical system associated with the equilibrium conditions.
    ChemicalSystem m_system;

    /// The names of the input variables in the chemical equilibrium calculation.
    Strings m_inputs;

    /// The *q* control variables in the chemical equilibrium problem.
    Vec<ControlVariableQ> qvars;

    /// The *p* control variables in the chemical equilibrium problem.
    Vec<ControlVariableP> pvars;

    /// The single equation constraints to be satisfied at chemical equilibrium.
    Vec<EquationConstraint> econstraints_single;

    /// The systems of equation constraints to be satisfied at chemical equilibrium.
    Vec<EquationConstraints> econstraints_system;

    /// The collected ids of every equation constraint above (single and system ones).
    Strings econstraints_ids;

    /// The single reactivity constraints to be satisfied at chemical equilibrium.
    Vec<ReactivityConstraint> rconstraints_single;

    /// The systems of reactivity constraints to be satisfied at chemical equilibrium.
    Vec<ReactivityConstraints> rconstraints_system;

    /// The collected ids of every reactivity constraint above (single and system ones).
    Strings rconstraints_ids;

    // ----- AUXILIARY DATA MEMBERS ----- //

    /// The chemical formulas of the explicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_explicit;

    /// The chemical formulas of the implicit titrants whose amounts are unknown.
    Vec<ChemicalFormula> titrants_implicit;

    /// The names of the species whose chemical potentials are unknown.
    Strings species_with_unknown_chemical_potentials;

private:
    /// Throw an error if a given titrant has already been registered explicitly or implicitly.
    /// @param substance The chemical formula of the titrant.
    auto throwErrorIfTitrantHasBeenRegistered(ChemicalFormula const& substance) const -> void;
};

} // namespace Reaktoro
