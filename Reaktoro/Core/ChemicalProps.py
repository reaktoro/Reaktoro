# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2024 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *
import pytest
from math import *


@pytest.fixture
def database() -> Database:
    return Database([
        Species("H2O(g)").withStandardGibbsEnergy(0.0),
        Species("CO2(g)").withStandardGibbsEnergy(0.0),
    ])


def testChemicalProps(database: Database) -> None:

    phases = Phases(database)
    phases.add( GaseousPhase("H2O(g) CO2(g)") )

    system = ChemicalSystem(phases)
    state = ChemicalState(system)

    state.setTemperature(100.0, "celsius")
    state.setPressure(1.0, "MPa")
    state.setSpeciesAmounts([3.0, 7.0])

    props = ChemicalProps(state)

    assert props.temperature() == state.temperature()
    assert props.pressure() == state.pressure()
    assert props.speciesAmounts() == state.speciesAmounts()
    assert props.speciesAmounts() == [3.0, 7.0]

    assert props.speciesMoleFractions() == [0.3, 0.7]
    assert props.speciesActivitiesLn().asarray() == pytest.approx([log(3.0), log(7.0)])
    # Partial molar volumes are only evaluated for phases with activity models with cubic EoSs.
    assert props.speciesPartialMolarVolumes() == 0.0


def testChemicalPropsPengRobinsonPartialMolarVolumes(database: Database) -> None:

    phases = Phases(database)
    phases.add( GaseousPhase("H2O(g) CO2(g)").set(ActivityModelPengRobinson()) )

    system = ChemicalSystem(phases)
    state = ChemicalState(system)

    state.setTemperature(0.0, "celsius")
    state.setPressure(1.0, "atm")
    state.setSpeciesAmounts([0.3, 0.7])

    props = ChemicalProps(state)

    assert props.temperature() == state.temperature()
    assert props.pressure() == state.pressure()
    assert props.speciesAmounts() == state.speciesAmounts()
    assert props.speciesAmounts() == [0.3, 0.7]

    assert props.speciesMoleFractions() == [0.3, 0.7]
    assert props.speciesPartialMolarVolumes().asarray() == pytest.approx([0.0220046, 0.0222578])
