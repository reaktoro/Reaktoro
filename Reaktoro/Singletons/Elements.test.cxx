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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Singletons/Elements.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Elements", "[Elements]")
{
    REQUIRE(Elements::size() == Elements::data().size());

    REQUIRE_NOTHROW( Elements::withSymbol("H"  ).value().name() == "Hydrogen"      );
    REQUIRE_NOTHROW( Elements::withSymbol("He" ).value().name() == "Helium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Li" ).value().name() == "Lithium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Be" ).value().name() == "Beryllium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("B"  ).value().name() == "Boron"         );
    REQUIRE_NOTHROW( Elements::withSymbol("C"  ).value().name() == "Carbon"        );
    REQUIRE_NOTHROW( Elements::withSymbol("N"  ).value().name() == "Nitrogen"      );
    REQUIRE_NOTHROW( Elements::withSymbol("O"  ).value().name() == "Oxygen"        );
    REQUIRE_NOTHROW( Elements::withSymbol("F"  ).value().name() == "Fluorine"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Ne" ).value().name() == "Neon"          );
    REQUIRE_NOTHROW( Elements::withSymbol("Na" ).value().name() == "Sodium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Mg" ).value().name() == "Magnesium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Al" ).value().name() == "Aluminum"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Si" ).value().name() == "Silicon"       );
    REQUIRE_NOTHROW( Elements::withSymbol("P"  ).value().name() == "Phosphorus"    );
    REQUIRE_NOTHROW( Elements::withSymbol("S"  ).value().name() == "Sulfur"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Cl" ).value().name() == "Chlorine"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Ar" ).value().name() == "Argon"         );
    REQUIRE_NOTHROW( Elements::withSymbol("K"  ).value().name() == "Potassium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Ca" ).value().name() == "Calcium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Sc" ).value().name() == "Scandium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Ti" ).value().name() == "Titanium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("V"  ).value().name() == "Vanadium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Cr" ).value().name() == "Chromium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Mn" ).value().name() == "Manganese"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Fe" ).value().name() == "Iron"          );
    REQUIRE_NOTHROW( Elements::withSymbol("Co" ).value().name() == "Cobalt"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Ni" ).value().name() == "Nickel"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Cu" ).value().name() == "Copper"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Zn" ).value().name() == "Zinc"          );
    REQUIRE_NOTHROW( Elements::withSymbol("Ga" ).value().name() == "Gallium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Ge" ).value().name() == "Germanium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("As" ).value().name() == "Arsenic"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Se" ).value().name() == "Selenium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Br" ).value().name() == "Bromine"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Kr" ).value().name() == "Krypton"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Rb" ).value().name() == "Rubidium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Sr" ).value().name() == "Strontium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Y"  ).value().name() == "Yttrium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Zr" ).value().name() == "Zirconium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Nb" ).value().name() == "Niobium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Mo" ).value().name() == "Molybdenum"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Tc" ).value().name() == "Technetium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Ru" ).value().name() == "Ruthenium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Rh" ).value().name() == "Rhodium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Pd" ).value().name() == "Palladium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Ag" ).value().name() == "Silver"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Cd" ).value().name() == "Cadmium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("In" ).value().name() == "Indium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Sn" ).value().name() == "Tin"           );
    REQUIRE_NOTHROW( Elements::withSymbol("Sb" ).value().name() == "Antimony"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Te" ).value().name() == "Tellurium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("I"  ).value().name() == "Iodine"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Xe" ).value().name() == "Xenon"         );
    REQUIRE_NOTHROW( Elements::withSymbol("Cs" ).value().name() == "Cesium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Ba" ).value().name() == "Barium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("La" ).value().name() == "Lanthanum"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Ce" ).value().name() == "Cerium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Pr" ).value().name() == "Praseodymium"  );
    REQUIRE_NOTHROW( Elements::withSymbol("Nd" ).value().name() == "Neodymium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Pm" ).value().name() == "Promethium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Sm" ).value().name() == "Samarium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Eu" ).value().name() == "Europium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Gd" ).value().name() == "Gadolinium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Tb" ).value().name() == "Terbium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Dy" ).value().name() == "Dysprosium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Ho" ).value().name() == "Holmium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Er" ).value().name() == "Erbium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Tm" ).value().name() == "Thulium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Yb" ).value().name() == "Ytterbium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Lu" ).value().name() == "Lutetium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Hf" ).value().name() == "Hafnium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Ta" ).value().name() == "Tantalum"      );
    REQUIRE_NOTHROW( Elements::withSymbol("W"  ).value().name() == "Tungsten"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Re" ).value().name() == "Rhenium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Os" ).value().name() == "Osmium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Ir" ).value().name() == "Iridium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Pt" ).value().name() == "Platinum"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Au" ).value().name() == "Gold"          );
    REQUIRE_NOTHROW( Elements::withSymbol("Hg" ).value().name() == "Mercury"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Tl" ).value().name() == "Thallium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Pb" ).value().name() == "Lead"          );
    REQUIRE_NOTHROW( Elements::withSymbol("Bi" ).value().name() == "Bismuth"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Po" ).value().name() == "Polonium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("At" ).value().name() == "Astatine"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Rn" ).value().name() == "Radon"         );
    REQUIRE_NOTHROW( Elements::withSymbol("Fr" ).value().name() == "Francium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Ra" ).value().name() == "Radium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Ac" ).value().name() == "Actinium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Th" ).value().name() == "Thorium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Pa" ).value().name() == "Protactinium"  );
    REQUIRE_NOTHROW( Elements::withSymbol("U"  ).value().name() == "Uranium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Np" ).value().name() == "Neptunium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Pu" ).value().name() == "Plutonium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Am" ).value().name() == "Americium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Cm" ).value().name() == "Curium"        );
    REQUIRE_NOTHROW( Elements::withSymbol("Bk" ).value().name() == "Berkelium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Cf" ).value().name() == "Californium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("Es" ).value().name() == "Einsteinium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("Fm" ).value().name() == "Fermium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Md" ).value().name() == "Mendelevium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("No" ).value().name() == "Nobelium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Lr" ).value().name() == "Lawrencium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Rf" ).value().name() == "Rutherfordium" );
    REQUIRE_NOTHROW( Elements::withSymbol("Db" ).value().name() == "Dubnium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Sg" ).value().name() == "Seaborgium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Bh" ).value().name() == "Bohrium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Hs" ).value().name() == "Hassium"       );
    REQUIRE_NOTHROW( Elements::withSymbol("Mt" ).value().name() == "Meitnerium"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Ds" ).value().name() == "Darmstadtium"  );
    REQUIRE_NOTHROW( Elements::withSymbol("Rg" ).value().name() == "Roentgenium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("Cn" ).value().name() == "Copernicium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("Nh" ).value().name() == "Nihonium"      );
    REQUIRE_NOTHROW( Elements::withSymbol("Fl" ).value().name() == "Flerovium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Mc" ).value().name() == "Moscovium"     );
    REQUIRE_NOTHROW( Elements::withSymbol("Lv" ).value().name() == "Livermorium"   );
    REQUIRE_NOTHROW( Elements::withSymbol("Ts" ).value().name() == "Tennessine"    );
    REQUIRE_NOTHROW( Elements::withSymbol("Og" ).value().name() == "Oganesson"     );

    REQUIRE_NOTHROW(Elements::withName("Hydrogen").value().symbol() == "H" );
    REQUIRE_NOTHROW(Elements::withName("Helium"  ).value().symbol() == "He");
    REQUIRE_NOTHROW(Elements::withName("Lithium" ).value().symbol() == "Li");

    REQUIRE_FALSE(Elements::withSymbol("Ab"));
    REQUIRE_FALSE(Elements::withSymbol("Xy"));

    Elements::append(Element().withSymbol("Aa").withTags({"tag1", "tag2", "tag3"}));
    Elements::append(Element().withSymbol("Bb").withTags({"tag1", "tag3"}));
    Elements::append(Element().withSymbol("Cc").withTags({"tag1"}));

    const auto elements_with_tag = Elements::withTag("tag1");

    REQUIRE(elements_with_tag.size() == 3);
    REQUIRE(elements_with_tag[0].symbol() == "Aa");
    REQUIRE(elements_with_tag[1].symbol() == "Bb");
    REQUIRE(elements_with_tag[2].symbol() == "Cc");

    const auto elements_with_tags = Elements::withTags({"tag3", "tag1"});

    REQUIRE(elements_with_tags.size() == 2);
    REQUIRE(elements_with_tags[0].symbol() == "Aa");
    REQUIRE(elements_with_tags[1].symbol() == "Bb");

    Elements::replace("Na", Element().withSymbol("New"));

    REQUIRE(Elements::withSymbol("Na").has_value() == false);
    REQUIRE(Elements::withSymbol("New").has_value());
    REQUIRE(Elements::withSymbol("New").value().symbol() == "New");

    Elements::clear();

    REQUIRE(Elements::size() == 0);

    Elements::reset();

    REQUIRE(Elements::withSymbol("Na"));
}
