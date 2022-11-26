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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

#include <reaktplot/reaktplot.hpp>
using namespace reaktplot;

int main()
{
    ReactionRateModelParamsPalandriKharaka params;
    params.mineral = "Calcite";
    params.mechanisms = {{
        { "Acid", -0.30, 14.4, 1.0, 1.0, {{"H+", "a", 1.0}} },
        { "Neutral", -5.81, 23.5, 1.0, 1.0 },
        { "Carbonate", -3.48, 35.4, 1.0, 1.0, {{"CO2", "P", 1.0}} }
    }};

    PhreeqcDatabase db("phreeqc.dat");

    ChemicalSystem system(db,
        AqueousPhase("H2O(aq) H+ OH- Ca+2 HCO3- CO3-2 CO2(aq)").set(ActivityModelDavies()),
        MineralPhase("Calcite"),
        MineralReaction("Calcite").setRateModel(ReactionRateModelPalandriKharaka(params)),
        MineralSurface("Calcite", 5.0, "cm2", 70, "mg", 0.667)  // surface area = (initial surface area)*((current mass)/(initial mass))**power = (5 cm2)*((current mass in mg)/(70 mg))**0.667
    );

    ChemicalState state(system);
    state.set("H2O(aq)", 1.0, "kg");
    state.set("Calcite", 70, "mg");

    AqueousProps aprops(system);

    KineticsSolver solver(system);

    Table table;

    auto dt = 2.0;  // time step (in seconds)

    for(auto i = 0; i < 501; ++i)
    {
        auto result = solver.solve(state, dt);

        errorif(result.failed(), "Calculation did not succeed at time step #", i, ".");

        aprops.update(state.props());

        table.column("Timestep") << i;
        table.column("Time")     << i*dt / 60.0; // from seconds to minute
        table.column("Calcite")  << state.props().speciesMass("Calcite") * 1e6; // from kg to mg
        table.column("Ca+2")     << aprops.speciesMolality("Ca+2");
        table.column("HCO3-")    << aprops.speciesMolality("HCO3-");
        table.column("CO3-2")    << aprops.speciesMolality("CO3-2");
        table.column("CO2(aq)")  << aprops.speciesMolality("CO2(aq)");
        table.column("pH")       << aprops.pH();
    }

    table.save("ex-kinetics-calcite-using-embedded-params.txt");

    Figure fig1;
    fig1.title("AQUEOUS SPECIES CONCENTRATIONS OVER TIME");
    fig1.xaxisTitle("Time [minute]");
    fig1.yaxisTitle("Concentration [molal]");
    fig1.drawLine(table["Time"], table["Ca+2"], "Ca<sup>+2</sup>");
    fig1.drawLine(table["Time"], table["HCO3-"], "HCO<sub>3</sub><sup>-</sup>");
    fig1.drawLine(table["Time"], table["CO3-2"], "CO<sub>3</sub><sup>-2</sup>");
    fig1.drawLine(table["Time"], table["CO2(aq)"], "CO<sub>2</sub>(aq)");
    fig1.yaxisScaleLog();
    fig1.show();
    fig1.save("ex-kinetics-calcite-using-embedded-params-fig1.pdf");

    Figure fig2;
    fig2.title("CALCITE MASS OVER TIME");
    fig2.xaxisTitle("Time [minute]");
    fig2.yaxisTitle("Mass [mg]");
    fig2.drawLine(table["Time"], table["Calcite"], "Calcite");
    fig2.show();
    fig2.save("ex-kinetics-calcite-using-embedded-params-fig2.pdf");

    Figure fig3;
    fig3.title("PH OVER TIME");
    fig3.xaxisTitle("Time [minute]");
    fig3.yaxisTitle("pH");
    fig3.drawLine(table["Time"], table["pH"], "pH");
    fig3.show();
    fig3.save("ex-kinetics-calcite-using-embedded-params-fig3.pdf");
}
