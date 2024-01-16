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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// The class that computes the full pivoting Auxiliary struct for storing the LU decomposition of a matrix `A`.
struct LU
{
    /// Construct a default LU instance.
    LU();

    /// Construct a LU instance with given matrix.
    explicit LU(MatrixXdConstRef A);

    /// Construct a LU instance with given matrix and scaling column-weights.
    LU(MatrixXdConstRef A, VectorXdConstRef W);

    /// Return true if empty.
    auto empty() const -> bool;

    /// Compute the LU decomposition of the given matrix.
    auto compute(MatrixXdConstRef A) -> void;

    /// Compute the LU decomposition of the given matrix with scaling column-weights.
    auto compute(MatrixXdConstRef A, VectorXdConstRef W) -> void;

    /// Solve the linear system `AX = B` using the calculated LU decomposition.
    auto solve(MatrixXdConstRef b) -> MatrixXd;

    /// Solve the linear system `tr(A)X = B` using the calculated LU decomposition.
    auto trsolve(MatrixXdConstRef B) -> MatrixXd;

    /// The last decomposed matrix A
    MatrixXd A_last;

    /// The last weights used for column scaling
    VectorXd W_last;

    /// The lower triangular matrix `L` in the LU decomposition of the matrix `PAQ = LU`.
    MatrixXd L;

    /// The upper triangular matrix `U` in the LU decomposition of the matrix `PAQ = LU`.
    MatrixXd U;

    /// The permutation matrix `P` in the LU decomposition of the matrix `PAQ = LU`.
    PermutationMatrix P;

    /// The permutation matrix `Q` in the LU decomposition of the matrix `PAQ = LU`.
    PermutationMatrix Q;

    /// The rank of the matrix `A`
    Index rank;
};

} // namespace Reaktoro
