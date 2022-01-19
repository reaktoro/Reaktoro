# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
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


def testChemicalProps():
    db = Database([
        Species("H2O(g)").withStandardGibbsEnergy(0.0),
        Species("CO2(g)").withStandardGibbsEnergy(0.0),
    ])

    phases = Phases(db)
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

