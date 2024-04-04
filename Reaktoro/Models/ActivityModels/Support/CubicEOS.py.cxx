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
#include <Reaktoro/Core/Model.py.hxx>
#include <Reaktoro/Models/ActivityModels/Support/CubicEOS.hpp>
using namespace Reaktoro;

void exportCubicEOS(py::module& m)
{
    auto ceos = m.def_submodule("CubicEOS");

    py::class_<CubicEOS::Substance>(ceos, "Substance")
        .def(py::init<String, real, real, real>(),
            py::arg("formula") = "",
            py::arg("Tcr") = NaN,
            py::arg("Pcr") = NaN,
            py::arg("omega") = NaN)
        .def_readwrite("formula", &CubicEOS::Substance::formula, "The chemical formula of the substance in the fluid phase.")
        .def_readwrite("Tcr", &CubicEOS::Substance::Tcr, "The critical temperature of the substance in the fluid phase (in K).")
        .def_readwrite("Pcr", &CubicEOS::Substance::Pcr, "The critical pressure of the substance in the fluid phase (in Pa).")
        .def_readwrite("omega", &CubicEOS::Substance::omega, "The acentric factor of the substance in the fluid phase.")
        ;

    py::class_<CubicEOS::Props>(ceos, "Props")
        .def(py::init<>())
        .def_readwrite("V", &CubicEOS::Props::V, "The molar volume of the phase (in m3/mol).")
        .def_readwrite("VT", &CubicEOS::Props::VT, "The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).")
        .def_readwrite("VP", &CubicEOS::Props::VP, "The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).")
        .def_readwrite("Gres", &CubicEOS::Props::Gres, "The residual molar Gibbs energy of the phase (in J/mol).")
        .def_readwrite("Hres", &CubicEOS::Props::Hres, "The residual molar enthalpy of the phase (in J/mol).")
        .def_readwrite("Cpres", &CubicEOS::Props::Cpres, "The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).")
        .def_readwrite("Cvres", &CubicEOS::Props::Cvres, "The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).")
        .def_readwrite("ln_phi", &CubicEOS::Props::ln_phi, "The ln fugacity coefficients of the species in the phase.")
        .def_readwrite("som", &CubicEOS::Props::som, "The state of matter of the fluid phase")
        ;

    auto createBipFromSize = [](Index size)
    {
        auto zero = MatrixXr::Zero(size, size);
        return CubicEOS::Bip{zero, zero, zero};
    };

    py::class_<CubicEOS::Bip>(ceos, "Bip")
        .def(py::init<>())
        .def(py::init(createBipFromSize))
        .def(py::init<MatrixXr, MatrixXr, MatrixXr>())
        .def_readwrite("k", &CubicEOS::Bip::k, "The computed binary interaction parameters kij for the species in the fluid phase.")
        .def_readwrite("kT", &CubicEOS::Bip::kT, "The computed first-order temperature derivative of the binary interaction parameters kij.")
        .def_readwrite("kTT", &CubicEOS::Bip::kTT, "The computed second-order temperature derivative of the binary interaction parameters kij.")
        ;

    py::class_<CubicEOS::BipModelArgs>(ceos, "BipModelArgs")
        .def_property_readonly("substances", [](CubicEOS::BipModelArgs const& self) { return self.substances; }, "The chemical formulas of the substances in the fluid phase")
        .def_property_readonly("T"         , [](CubicEOS::BipModelArgs const& self) { return self.T; }         , "The temperature of the fluid (in K)")
        .def_property_readonly("Tcr"       , [](CubicEOS::BipModelArgs const& self) { return self.Tcr; }       , "The critical temperatures of the substances in the fluid phase (in K)")
        .def_property_readonly("Pcr"       , [](CubicEOS::BipModelArgs const& self) { return self.Pcr; }       , "The critical pressures of the substances in the fluid phase (in Pa)")
        .def_property_readonly("omega"     , [](CubicEOS::BipModelArgs const& self) { return self.omega; }     , "The acentric factors of the substances in the fluid phase")
        .def_property_readonly("a"         , [](CubicEOS::BipModelArgs const& self) { return self.a; }         , "The variables a[i] in the cubic equation of state for each substance in the fluid phase computed at current temperature")
        .def_property_readonly("aT"        , [](CubicEOS::BipModelArgs const& self) { return self.aT; }        , "The first-order temperature derivative of a[i]")
        .def_property_readonly("aTT"       , [](CubicEOS::BipModelArgs const& self) { return self.aTT; }       , "The second-order temperature derivative of a[i]")
        .def_property_readonly("alpha"     , [](CubicEOS::BipModelArgs const& self) { return self.alpha; }     , "The variables alpha[i] in the cubic equation of state for each substance in the fluid phase computed at current temperature")
        .def_property_readonly("alphaT"    , [](CubicEOS::BipModelArgs const& self) { return self.alphaT; }    , "The first-order temperature derivative of alpha[i]")
        .def_property_readonly("alphaTT"   , [](CubicEOS::BipModelArgs const& self) { return self.alphaTT; }   , "The second-order temperature derivative of alpha[i]")
        .def_property_readonly("b"         , [](CubicEOS::BipModelArgs const& self) { return self.b; }         , "The variables b[i] in the cubic equation of state for each substance in the fluid phase computed at current temperature")
        ;

    exportModel<CubicEOS::Bip, CubicEOS::BipModelArgs const&>(ceos, "BipModel");

    py::class_<CubicEOS::Alpha>(ceos, "Alpha")
        .def(py::init<real, real, real>(),
            py::arg("alpha") = 0.0,
            py::arg("alphaT") = 0.0,
            py::arg("alphaTT") = 0.0)
        .def_readwrite("alpha", &CubicEOS::Alpha::alpha, "The evaluated value of the alpha(T; w) function in a cubic equation of state.")
        .def_readwrite("alphaT", &CubicEOS::Alpha::alphaT, "The evaluated first-order temperature derivative of the alpha(T; w) function.")
        .def_readwrite("alphaTT", &CubicEOS::Alpha::alphaTT, "The evaluated second-order temperature derivative of the alpha(T; w) function.")
        ;

    py::class_<CubicEOS::AlphaModelArgs>(ceos, "AlphaModelArgs")
        .def_property_readonly("Tr", [](CubicEOS::AlphaModelArgs const& self) { return self.Tr; }, "The reduced temperature Tr = T/Tcr with respect to the substance for which alpha(T; w) is computed.")
        .def_property_readonly("TrT", [](CubicEOS::AlphaModelArgs const& self) { return self.TrT; }, "The first-order temperature derivative of Tr for convenience.")
        .def_property_readonly("omega", [](CubicEOS::AlphaModelArgs const& self) { return self.omega; }, "The acentric factor of the substance for which alpha(T; w) is computed.")
        ;

    exportModel<CubicEOS::Alpha, CubicEOS::AlphaModelArgs const&>(ceos, "AlphaModel");

    ceos.def("EquationModelVanDerWaals", &CubicEOS::EquationModelVanDerWaals, "Return a cubic equation model representative of the van der Waals (1873) cubic equation of state.");
    ceos.def("EquationModelRedlichKwong", &CubicEOS::EquationModelRedlichKwong, "Return a cubic equation model representative of the Redlich-Kwong (1949) cubic equation of state.");
    ceos.def("EquationModelSoaveRedlichKwong", &CubicEOS::EquationModelSoaveRedlichKwong, "Return a cubic equation model representative of the Soave-Redlich-Kwong (1972) cubic equation of state.");
    ceos.def("EquationModelPengRobinson76", &CubicEOS::EquationModelPengRobinson76, "Return a cubic equation model representative of the Peng-Robinson (1976) cubic equation of state.");
    ceos.def("EquationModelPengRobinson78", &CubicEOS::EquationModelPengRobinson78, "Return a cubic equation model representative of the Peng-Robinson (1978) cubic equation of state.");
    ceos.def("EquationModelPengRobinson", &CubicEOS::EquationModelPengRobinson, "Return a cubic equation model representative of the Peng-Robinson (1978) cubic equation of state.");

    py::class_<CubicEOS::EquationSpecs>(ceos, "EquationSpecs")
        .def_readwrite("substances", &CubicEOS::EquationSpecs::substances, "The substances in the fluid phase and their attributes.")
        .def_readwrite("eqmodel", &CubicEOS::EquationSpecs::eqmodel, "The cubic equation of state model to be used.")
        .def_readwrite("bipmodel", &CubicEOS::EquationSpecs::bipmodel, "The function that calculates binary interaction parameters kij.")
        ;

    py::class_<CubicEOS::Equation>(ceos, "Equation")
        .def(py::init<CubicEOS::EquationSpecs>())
        .def("equationSpecs", &CubicEOS::Equation::equationSpecs, "Return the underlying EquationSpecs object used to create this Equation object.")
        .def("compute", &CubicEOS::Equation::compute, "Compute the thermodynamic properties of the phase.")
        ;

    py::class_<CubicEOS::BipModelParamsPhreeqc>(ceos, "BipModelParamsPhreeqc")
        .def_readwrite("kH2O_CO2", &CubicEOS::BipModelParamsPhreeqc::kH2O_CO2, "The binary interaction parameter k_ij for the substance pair H2O-CO2.")
        .def_readwrite("kH2O_H2S", &CubicEOS::BipModelParamsPhreeqc::kH2O_H2S, "The binary interaction parameter k_ij for the substance pair H2O-H2S.")
        .def_readwrite("kH2O_CH4", &CubicEOS::BipModelParamsPhreeqc::kH2O_CH4, "The binary interaction parameter k_ij for the substance pair H2O-CH4.")
        .def_readwrite("kH2O_N2", &CubicEOS::BipModelParamsPhreeqc::kH2O_N2, "The binary interaction parameter k_ij for the substance pair H2O-N2.")
        ;

    py::class_<CubicEOS::BipModelParamsSoreideWhitson>(ceos, "BipModelParamsSoreideWhitson")
        .def_readwrite("kH2O_CO2", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_CO2, "The binary interaction parameter k_ij for the substance pair H2O/CO2.")
        .def_readwrite("kH2O_N2", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_N2, "The binary interaction parameter k_ij for the substance pair H2O/N2.")
        .def_readwrite("kH2O_CH4", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_CH4, "The binary interaction parameter k_ij for the substance pair H2O/CH4.")
        .def_readwrite("kH2O_C2H6", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_C2H6, "The binary interaction parameter k_ij for the substance pair H2O/C2H6.")
        .def_readwrite("kH2O_C3H8", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_C3H8, "The binary interaction parameter k_ij for the substance pair H2O/C3H8.")
        .def_readwrite("kH2O_nC4H10", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_nC4H10, "The binary interaction parameter k_ij for the substance pair H2O/n-C4H10.")
        .def_readwrite("kH2O_H2S_a1", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_H2S_a1, "The coefficient a_1 when computing the binary interaction parameter k_ij for the substance pair H2O/H2S.")
        .def_readwrite("kH2O_H2S_a2", &CubicEOS::BipModelParamsSoreideWhitson::kH2O_H2S_a2, "The coefficient a_2 when computing the binary interaction parameter k_ij for the substance pair H2O/H2S.")
        ;

    ceos.def("BipModelPhreeqc", &CubicEOS::BipModelPhreeqc, "Return a binary interaction parameter model for Peng-Robinson EOS equivalent to that used in PHREEQC.");
    ceos.def("BipModelSoreideWhitson", &CubicEOS::BipModelSoreideWhitson, "Return a binary interaction parameter model for Peng-Robinson EOS equivalent to that reported in Søreide and Whitson (1992).");

    // DEPRECATED METHODS TO BE REMOVED IN THE NEAR FUTURE

    struct DeprecatedCubicEOSBipModelParamsPHREEQC {};
    py::class_<DeprecatedCubicEOSBipModelParamsPHREEQC>(ceos, "BipModelParamsPHREEQC")
        .def(py::init([]() -> DeprecatedCubicEOSBipModelParamsPHREEQC { errorif(true, "BipModelPHREEQC has been deprecated and renamed to BipModelPhreeqc."); return {}; } ))
        ;

    ceos.def("BipModelPHREEQC", [](Strings const&, CubicEOS::BipModelParamsPhreeqc const&) -> CubicEOS::BipModel { errorif(true, "BipModelPHREEQC has been deprecated and renamed to BipModelPhreeqc."); return {}; },
        "BipModelPHREEQC has been deprecated and renamed to BipModelPhreeqc.");
}
