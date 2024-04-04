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

#include "Elements.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace detail {

/// The default elements for the Elements object.
const Vec<Element> default_elements =
{
    Element({ "H"  , 0.001007940 , "Hydrogen"      }),
    Element({ "He" , 0.004002602 , "Helium"        }),
    Element({ "Li" , 0.006941000 , "Lithium"       }),
    Element({ "Be" , 0.009012180 , "Beryllium"     }),
    Element({ "B"  , 0.010811000 , "Boron"         }),
    Element({ "C"  , 0.012011000 , "Carbon"        }),
    Element({ "N"  , 0.014006740 , "Nitrogen"      }),
    Element({ "O"  , 0.015999400 , "Oxygen"        }),
    Element({ "F"  , 0.018998403 , "Fluorine"      }),
    Element({ "Ne" , 0.020179700 , "Neon"          }),
    Element({ "Na" , 0.022989768 , "Sodium"        }),
    Element({ "Mg" , 0.024305000 , "Magnesium"     }),
    Element({ "Al" , 0.026981539 , "Aluminum"      }),
    Element({ "Si" , 0.028085500 , "Silicon"       }),
    Element({ "P"  , 0.030973762 , "Phosphorus"    }),
    Element({ "S"  , 0.032066000 , "Sulfur"        }),
    Element({ "Cl" , 0.035452700 , "Chlorine"      }),
    Element({ "Ar" , 0.039948000 , "Argon"         }),
    Element({ "K"  , 0.039098300 , "Potassium"     }),
    Element({ "Ca" , 0.040078000 , "Calcium"       }),
    Element({ "Sc" , 0.044955910 , "Scandium"      }),
    Element({ "Ti" , 0.047880000 , "Titanium"      }),
    Element({ "V"  , 0.050941500 , "Vanadium"      }),
    Element({ "Cr" , 0.051996100 , "Chromium"      }),
    Element({ "Mn" , 0.054938050 , "Manganese"     }),
    Element({ "Fe" , 0.055847000 , "Iron"          }),
    Element({ "Co" , 0.058933200 , "Cobalt"        }),
    Element({ "Ni" , 0.058693400 , "Nickel"        }),
    Element({ "Cu" , 0.063546000 , "Copper"        }),
    Element({ "Zn" , 0.065390000 , "Zinc"          }),
    Element({ "Ga" , 0.069723000 , "Gallium"       }),
    Element({ "Ge" , 0.072610000 , "Germanium"     }),
    Element({ "As" , 0.074921590 , "Arsenic"       }),
    Element({ "Se" , 0.078960000 , "Selenium"      }),
    Element({ "Br" , 0.079904000 , "Bromine"       }),
    Element({ "Kr" , 0.083800000 , "Krypton"       }),
    Element({ "Rb" , 0.085467800 , "Rubidium"      }),
    Element({ "Sr" , 0.087620000 , "Strontium"     }),
    Element({ "Y"  , 0.088905850 , "Yttrium"       }),
    Element({ "Zr" , 0.091224000 , "Zirconium"     }),
    Element({ "Nb" , 0.092906380 , "Niobium"       }),
    Element({ "Mo" , 0.095940000 , "Molybdenum"    }),
    Element({ "Tc" , 0.097907200 , "Technetium"    }),
    Element({ "Ru" , 0.101070000 , "Ruthenium"     }),
    Element({ "Rh" , 0.102905500 , "Rhodium"       }),
    Element({ "Pd" , 0.106420000 , "Palladium"     }),
    Element({ "Ag" , 0.107868200 , "Silver"        }),
    Element({ "Cd" , 0.112411000 , "Cadmium"       }),
    Element({ "In" , 0.114818000 , "Indium"        }),
    Element({ "Sn" , 0.118710000 , "Tin"           }),
    Element({ "Sb" , 0.121760000 , "Antimony"      }),
    Element({ "Te" , 0.127600000 , "Tellurium"     }),
    Element({ "I"  , 0.126904470 , "Iodine"        }),
    Element({ "Xe" , 0.131290000 , "Xenon"         }),
    Element({ "Cs" , 0.132905430 , "Cesium"        }),
    Element({ "Ba" , 0.137327000 , "Barium"        }),
    Element({ "La" , 0.138905500 , "Lanthanum"     }),
    Element({ "Ce" , 0.140115000 , "Cerium"        }),
    Element({ "Pr" , 0.140907650 , "Praseodymium"  }),
    Element({ "Nd" , 0.144240000 , "Neodymium"     }),
    Element({ "Pm" , 0.144912700 , "Promethium"    }),
    Element({ "Sm" , 0.150360000 , "Samarium"      }),
    Element({ "Eu" , 0.151965000 , "Europium"      }),
    Element({ "Gd" , 0.157250000 , "Gadolinium"    }),
    Element({ "Tb" , 0.158925340 , "Terbium"       }),
    Element({ "Dy" , 0.162500000 , "Dysprosium"    }),
    Element({ "Ho" , 0.164930320 , "Holmium"       }),
    Element({ "Er" , 0.167260000 , "Erbium"        }),
    Element({ "Tm" , 0.168934210 , "Thulium"       }),
    Element({ "Yb" , 0.173040000 , "Ytterbium"     }),
    Element({ "Lu" , 0.174967000 , "Lutetium"      }),
    Element({ "Hf" , 0.178490000 , "Hafnium"       }),
    Element({ "Ta" , 0.180947900 , "Tantalum"      }),
    Element({ "W"  , 0.183840000 , "Tungsten"      }),
    Element({ "Re" , 0.186207000 , "Rhenium"       }),
    Element({ "Os" , 0.190230000 , "Osmium"        }),
    Element({ "Ir" , 0.192220000 , "Iridium"       }),
    Element({ "Pt" , 0.195080000 , "Platinum"      }),
    Element({ "Au" , 0.196966540 , "Gold"          }),
    Element({ "Hg" , 0.200590000 , "Mercury"       }),
    Element({ "Tl" , 0.204383300 , "Thallium"      }),
    Element({ "Pb" , 0.207200000 , "Lead"          }),
    Element({ "Bi" , 0.208980370 , "Bismuth"       }),
    Element({ "Po" , 0.208982400 , "Polonium"      }),
    Element({ "At" , 0.209987100 , "Astatine"      }),
    Element({ "Rn" , 0.222017600 , "Radon"         }),
    Element({ "Fr" , 0.223019700 , "Francium"      }),
    Element({ "Ra" , 0.226025400 , "Radium"        }),
    Element({ "Ac" , 0.227027800 , "Actinium"      }),
    Element({ "Th" , 0.232038100 , "Thorium"       }),
    Element({ "Pa" , 0.231035880 , "Protactinium"  }),
    Element({ "U"  , 0.238028900 , "Uranium"       }),
    Element({ "Np" , 0.237048000 , "Neptunium"     }),
    Element({ "Pu" , 0.244064200 , "Plutonium"     }),
    Element({ "Am" , 0.243061400 , "Americium"     }),
    Element({ "Cm" , 0.247070300 , "Curium"        }),
    Element({ "Bk" , 0.247070300 , "Berkelium"     }),
    Element({ "Cf" , 0.251079600 , "Californium"   }),
    Element({ "Es" , 0.252083000 , "Einsteinium"   }),
    Element({ "Fm" , 0.257095100 , "Fermium"       }),
    Element({ "Md" , 0.258100000 , "Mendelevium"   }),
    Element({ "No" , 0.259100900 , "Nobelium"      }),
    Element({ "Lr" , 0.262110000 , "Lawrencium"    }),
    Element({ "Rf" , 0.261000000 , "Rutherfordium" }),
    Element({ "Db" , 0.262000000 , "Dubnium"       }),
    Element({ "Sg" , 0.266000000 , "Seaborgium"    }),
    Element({ "Bh" , 0.264000000 , "Bohrium"       }),
    Element({ "Hs" , 0.269000000 , "Hassium"       }),
    Element({ "Mt" , 0.268000000 , "Meitnerium"    }),
    Element({ "Ds" , 0.269000000 , "Darmstadtium"  }),
    Element({ "Rg" , 0.272000000 , "Roentgenium"   }),
    Element({ "Cn" , 0.277000000 , "Copernicium"   }),
    Element({ "Nh" , 0.000000000 , "Nihonium"      }),
    Element({ "Fl" , 0.289000000 , "Flerovium"     }),
    Element({ "Mc" , 0.000000000 , "Moscovium"     }),
    Element({ "Lv" , 0.000000000 , "Livermorium"   }),
    Element({ "Ts" , 0.000000000 , "Tennessine"    }),
    Element({ "Og" , 0.000000000 , "Oganesson"     }),
    Element({ "D"  , 0.002014102 , "Deuterium"     }),
    Element({ "T"  , 0.003016049 , "Tritium"       }),
};

} // namespace detail

Elements::Elements()
: m_elements(detail::default_elements)
{}

Elements::~Elements()
{}

auto Elements::instance() -> Elements&
{
    static Elements obj;
    return obj;
}

auto Elements::data() -> Vec<Element> const&
{
    return instance().m_elements;
}

auto Elements::elements() -> Vec<Element> const&
{
    return instance().m_elements;
}

auto Elements::clear() -> void
{
    auto& elements = instance().m_elements;
    elements.clear();
}

auto Elements::reset() -> void
{
    auto& elements = instance().m_elements;
    elements = detail::default_elements;
}

auto Elements::append(Element element) -> void
{
    auto& elements = instance().m_elements;
    elements.emplace_back(std::move(element));
}

auto Elements::replace(String const& existingSymbol, Element const& replacement) -> bool
{
    auto& elements = instance().m_elements;
    auto it = std::find_if(elements.begin(), elements.end(),
        [&](auto&& e) { return e.symbol() == existingSymbol; });
    if(it == elements.end())
        return false;
    *it = std::move(replacement);
    return true;
}

auto Elements::size() -> std::size_t
{
    return data().size();
}

auto Elements::withSymbol(String symbol) -> Optional<Element>
{
    const auto idx = indexfn(data(), [&](auto&& e) { return e.symbol() == symbol; });
    if(idx < size()) return data()[idx];
    return {};
}

auto Elements::withName(String name) -> Optional<Element>
{
    const auto idx = indexfn(data(), [&](auto&& e) { return e.name() == name; });
    if(idx < size()) return data()[idx];
    return {};
}

auto Elements::withTag(String tag) -> Vec<Element>
{
    return filter(data(), [&](auto&& e) { return contains(e.tags(), tag); });
}

auto Elements::withTags(StringList const& tags) -> Vec<Element>
{
    return filter(data(), [&](auto&& e) { return contained(tags, e.tags()); });
}

} // namespace Reaktoro
