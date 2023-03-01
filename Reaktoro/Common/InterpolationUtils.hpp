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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Vec<double>& scalars) -> Fn<real(real, real)>;

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Fn<double(double, double)>& func) -> Fn<real(real, real)>;

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Vec<Fn<double(double, double)>>& fs) -> Fn<ArrayXr(real, real)>;

/// Calculate a linear interpolation of *y* at *x* with given pairs *(x0, y0)* and *(x1, y1)*.
template<typename T, typename X, typename Y>
auto interpolateLinear(T const& x, X const& x0, X const& x1, Y const& y0, Y const& y1)
{
    using R = std::decay_t<decltype(y0 + (y1 - y0)/(x1 - x0) * (x - x0))>;
    assert(x0 <= x1);
    if(x0 == x1) return static_cast<R>(y0);
    return y0 + (y1 - y0)/(x1 - x0) * (x - x0);
}

/// Calculate a quadratic interpolation of *y* at *x* with given pairs *(x0, y0)* *(x1, y1)* and *(x2, y2)*.
template<typename T, typename X, typename Y>
auto interpolateQuadratic(T const& x, X const& x0, X const& x1, X const& x2, Y const& y0, Y const& y1, Y const& y2)
{
    assert(x0 <= x1 && x1 <= x2);
    if(x0 == x1 || x1 == x2) return interpolateLinear(x, x0, x2, y0, y2);
    const auto l0 = ((x - x1)*(x - x2))/((x0 - x1)*(x0 - x2));
    const auto l1 = ((x - x0)*(x - x2))/((x1 - x0)*(x1 - x2));
    const auto l2 = ((x - x0)*(x - x1))/((x2 - x0)*(x2 - x1));
    return y0*l0 + y1*l1 + y2*l2;
}

} // namespace Reaktoro
