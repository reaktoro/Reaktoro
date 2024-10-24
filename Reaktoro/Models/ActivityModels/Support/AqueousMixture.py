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

import numpy as np
import pytest

from reaktoro import *


def rhofn(T, P):
    wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter.Liquid)
    return wtp.D


def epsilonfn(T, P):
    wtp = waterThermoPropsWagnerPrussMemoized(T, P, StateOfMatter.Liquid)
    wep = waterElectroPropsJohnsonNorton(T, P, wtp)
    return wep.epsilon


def testAqueousMixture():
    species = SpeciesList(
        "H2O H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO3-- K+ CO2 HCl NaCl NaOH CaCl2 MgCl2 CaCO3 MgCO3"
    )

    T = autodiff.real(345.67)
    P = autodiff.real(123.4e5)

    n = np.random.rand(species.size())
    n[0] = 55.508

    x = n / n.sum()

    mixture = AqueousMixture(species)
    state = mixture.state(T, P, x)

    assert state.rho == pytest.approx(997.0470390177028)
    assert state.epsilon == pytest.approx(78.2451448082024)

    AqueousMixture.setDefaultWaterDensityFn(rhofn)
    AqueousMixture.setDefaultWaterDielectricConstantFn(epsilonfn)

    mixture = AqueousMixture(species)
    state = mixture.state(T, P, x)

    assert state.rho == pytest.approx(981.650989015948)
    assert state.epsilon == pytest.approx(63.37949243846789)

    AqueousMixture.resetDefaultWaterDensityFn()
    AqueousMixture.resetDefaultWaterDielectricConstantFn()

    mixture = AqueousMixture(species)
    state = mixture.state(T, P, x)

    assert state.rho == pytest.approx(997.0470390177028)
    assert state.epsilon == pytest.approx(78.2451448082024)
