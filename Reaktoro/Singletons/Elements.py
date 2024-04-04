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


def testElements():

    assert Elements.size() == len(Elements.data())

    assert Elements.withSymbol("H"  ).name() == "Hydrogen"
    assert Elements.withSymbol("He" ).name() == "Helium"
    assert Elements.withSymbol("Li" ).name() == "Lithium"
    assert Elements.withSymbol("Be" ).name() == "Beryllium"
    assert Elements.withSymbol("B"  ).name() == "Boron"
    assert Elements.withSymbol("C"  ).name() == "Carbon"
    assert Elements.withSymbol("N"  ).name() == "Nitrogen"
    assert Elements.withSymbol("O"  ).name() == "Oxygen"
    assert Elements.withSymbol("F"  ).name() == "Fluorine"
    assert Elements.withSymbol("Ne" ).name() == "Neon"
    assert Elements.withSymbol("Na" ).name() == "Sodium"
    assert Elements.withSymbol("Mg" ).name() == "Magnesium"
    assert Elements.withSymbol("Al" ).name() == "Aluminum"
    assert Elements.withSymbol("Si" ).name() == "Silicon"
    assert Elements.withSymbol("P"  ).name() == "Phosphorus"
    assert Elements.withSymbol("S"  ).name() == "Sulfur"
    assert Elements.withSymbol("Cl" ).name() == "Chlorine"
    assert Elements.withSymbol("Ar" ).name() == "Argon"
    assert Elements.withSymbol("K"  ).name() == "Potassium"
    assert Elements.withSymbol("Ca" ).name() == "Calcium"
    assert Elements.withSymbol("Sc" ).name() == "Scandium"
    assert Elements.withSymbol("Ti" ).name() == "Titanium"
    assert Elements.withSymbol("V"  ).name() == "Vanadium"
    assert Elements.withSymbol("Cr" ).name() == "Chromium"
    assert Elements.withSymbol("Mn" ).name() == "Manganese"
    assert Elements.withSymbol("Fe" ).name() == "Iron"
    assert Elements.withSymbol("Co" ).name() == "Cobalt"
    assert Elements.withSymbol("Ni" ).name() == "Nickel"
    assert Elements.withSymbol("Cu" ).name() == "Copper"
    assert Elements.withSymbol("Zn" ).name() == "Zinc"
    assert Elements.withSymbol("Ga" ).name() == "Gallium"
    assert Elements.withSymbol("Ge" ).name() == "Germanium"
    assert Elements.withSymbol("As" ).name() == "Arsenic"
    assert Elements.withSymbol("Se" ).name() == "Selenium"
    assert Elements.withSymbol("Br" ).name() == "Bromine"
    assert Elements.withSymbol("Kr" ).name() == "Krypton"
    assert Elements.withSymbol("Rb" ).name() == "Rubidium"
    assert Elements.withSymbol("Sr" ).name() == "Strontium"
    assert Elements.withSymbol("Y"  ).name() == "Yttrium"
    assert Elements.withSymbol("Zr" ).name() == "Zirconium"
    assert Elements.withSymbol("Nb" ).name() == "Niobium"
    assert Elements.withSymbol("Mo" ).name() == "Molybdenum"
    assert Elements.withSymbol("Tc" ).name() == "Technetium"
    assert Elements.withSymbol("Ru" ).name() == "Ruthenium"
    assert Elements.withSymbol("Rh" ).name() == "Rhodium"
    assert Elements.withSymbol("Pd" ).name() == "Palladium"
    assert Elements.withSymbol("Ag" ).name() == "Silver"
    assert Elements.withSymbol("Cd" ).name() == "Cadmium"
    assert Elements.withSymbol("In" ).name() == "Indium"
    assert Elements.withSymbol("Sn" ).name() == "Tin"
    assert Elements.withSymbol("Sb" ).name() == "Antimony"
    assert Elements.withSymbol("Te" ).name() == "Tellurium"
    assert Elements.withSymbol("I"  ).name() == "Iodine"
    assert Elements.withSymbol("Xe" ).name() == "Xenon"
    assert Elements.withSymbol("Cs" ).name() == "Cesium"
    assert Elements.withSymbol("Ba" ).name() == "Barium"
    assert Elements.withSymbol("La" ).name() == "Lanthanum"
    assert Elements.withSymbol("Ce" ).name() == "Cerium"
    assert Elements.withSymbol("Pr" ).name() == "Praseodymium"
    assert Elements.withSymbol("Nd" ).name() == "Neodymium"
    assert Elements.withSymbol("Pm" ).name() == "Promethium"
    assert Elements.withSymbol("Sm" ).name() == "Samarium"
    assert Elements.withSymbol("Eu" ).name() == "Europium"
    assert Elements.withSymbol("Gd" ).name() == "Gadolinium"
    assert Elements.withSymbol("Tb" ).name() == "Terbium"
    assert Elements.withSymbol("Dy" ).name() == "Dysprosium"
    assert Elements.withSymbol("Ho" ).name() == "Holmium"
    assert Elements.withSymbol("Er" ).name() == "Erbium"
    assert Elements.withSymbol("Tm" ).name() == "Thulium"
    assert Elements.withSymbol("Yb" ).name() == "Ytterbium"
    assert Elements.withSymbol("Lu" ).name() == "Lutetium"
    assert Elements.withSymbol("Hf" ).name() == "Hafnium"
    assert Elements.withSymbol("Ta" ).name() == "Tantalum"
    assert Elements.withSymbol("W"  ).name() == "Tungsten"
    assert Elements.withSymbol("Re" ).name() == "Rhenium"
    assert Elements.withSymbol("Os" ).name() == "Osmium"
    assert Elements.withSymbol("Ir" ).name() == "Iridium"
    assert Elements.withSymbol("Pt" ).name() == "Platinum"
    assert Elements.withSymbol("Au" ).name() == "Gold"
    assert Elements.withSymbol("Hg" ).name() == "Mercury"
    assert Elements.withSymbol("Tl" ).name() == "Thallium"
    assert Elements.withSymbol("Pb" ).name() == "Lead"
    assert Elements.withSymbol("Bi" ).name() == "Bismuth"
    assert Elements.withSymbol("Po" ).name() == "Polonium"
    assert Elements.withSymbol("At" ).name() == "Astatine"
    assert Elements.withSymbol("Rn" ).name() == "Radon"
    assert Elements.withSymbol("Fr" ).name() == "Francium"
    assert Elements.withSymbol("Ra" ).name() == "Radium"
    assert Elements.withSymbol("Ac" ).name() == "Actinium"
    assert Elements.withSymbol("Th" ).name() == "Thorium"
    assert Elements.withSymbol("Pa" ).name() == "Protactinium"
    assert Elements.withSymbol("U"  ).name() == "Uranium"
    assert Elements.withSymbol("Np" ).name() == "Neptunium"
    assert Elements.withSymbol("Pu" ).name() == "Plutonium"
    assert Elements.withSymbol("Am" ).name() == "Americium"
    assert Elements.withSymbol("Cm" ).name() == "Curium"
    assert Elements.withSymbol("Bk" ).name() == "Berkelium"
    assert Elements.withSymbol("Cf" ).name() == "Californium"
    assert Elements.withSymbol("Es" ).name() == "Einsteinium"
    assert Elements.withSymbol("Fm" ).name() == "Fermium"
    assert Elements.withSymbol("Md" ).name() == "Mendelevium"
    assert Elements.withSymbol("No" ).name() == "Nobelium"
    assert Elements.withSymbol("Lr" ).name() == "Lawrencium"
    assert Elements.withSymbol("Rf" ).name() == "Rutherfordium"
    assert Elements.withSymbol("Db" ).name() == "Dubnium"
    assert Elements.withSymbol("Sg" ).name() == "Seaborgium"
    assert Elements.withSymbol("Bh" ).name() == "Bohrium"
    assert Elements.withSymbol("Hs" ).name() == "Hassium"
    assert Elements.withSymbol("Mt" ).name() == "Meitnerium"
    assert Elements.withSymbol("Ds" ).name() == "Darmstadtium"
    assert Elements.withSymbol("Rg" ).name() == "Roentgenium"
    assert Elements.withSymbol("Cn" ).name() == "Copernicium"
    assert Elements.withSymbol("Nh" ).name() == "Nihonium"
    assert Elements.withSymbol("Fl" ).name() == "Flerovium"
    assert Elements.withSymbol("Mc" ).name() == "Moscovium"
    assert Elements.withSymbol("Lv" ).name() == "Livermorium"
    assert Elements.withSymbol("Ts" ).name() == "Tennessine"
    assert Elements.withSymbol("Og" ).name() == "Oganesson"

    assert Elements.withName("Hydrogen").symbol() == "H"
    assert Elements.withName("Helium"  ).symbol() == "He"
    assert Elements.withName("Lithium" ).symbol() == "Li"

    assert Elements.withSymbol("Ab") == None
    assert Elements.withSymbol("Xy") == None

    Elements.append(Element().withSymbol("Aa").withTags(["tag1", "tag2", "tag3"]))
    Elements.append(Element().withSymbol("Bb").withTags(["tag1", "tag3"]))
    Elements.append(Element().withSymbol("Cc").withTags(["tag1"]))

    elements_with_tag = Elements.withTag("tag1")

    assert len(elements_with_tag) == 3
    assert elements_with_tag[0].symbol() == "Aa"
    assert elements_with_tag[1].symbol() == "Bb"
    assert elements_with_tag[2].symbol() == "Cc"

    elements_with_tags = Elements.withTags(["tag3", "tag1"])

    assert len(elements_with_tags) == 2
    assert elements_with_tags[0].symbol() == "Aa"
    assert elements_with_tags[1].symbol() == "Bb"

    Elements.replace("Na", Element().withSymbol("New"))

    assert Elements.withSymbol("Na") == None
    assert Elements.withSymbol("New").symbol() == "New"

    Elements.clear()

    assert Elements.size() == 0

    Elements.reset()

    assert Elements.withSymbol("Na") != None
