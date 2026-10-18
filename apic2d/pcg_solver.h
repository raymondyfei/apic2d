// Copyright 2016 Raymond Yun Fei, Christopher Batty, Robert Bridson
//
// Licensed under the Apache License,
// Version 2.0(the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef APIC2D_PCG_SOLVER_H_
#define APIC2D_PCG_SOLVER_H_

// Implements PCG with Modified Incomplete Cholesky (0) preconditioner.
// PCGSolver<T> is the main class for setting up and solving a linear system.
// Note that this only handles symmetric positive (semi-)definite matrices,
// with guarantees made only for M-matrices (where off-diagonal entries are all
// non-positive, and row sums are non-negative).

#include <chrono>
#include <cmath>

#include "blas_wrapper.h"
#include "sparse_matrix.h"

//============================================================================
// A simple compressed sparse column data structure (with separate diagonal)
// for lower triangular matrices
namespace robertbridson {
template <class T>
struct SparseColumnLowerFactor {
  unsigned int n;
  std::vector<T> invdiag;              // reciprocals of diagonal elements
  std::vector<T> value;                // values below the diagonal, listed column by column
  std::vector<unsigned int> rowindex;  // a list of all row indices, for each column in turn
  std::vector<unsigned int> colstart;  // where each column begins in rowindex (plus an extra entry at
                                       // the end, of #nonzeros)
  std::vector<T> adiag;                // just used in factorization: minimum "safe" diagonal
                                       // entry allowed

  explicit SparseColumnLowerFactor(unsigned int n_ = 0) : n(n_), invdiag(n_), colstart(n_ + 1), adiag(n_) {}

  void clear(void) {
    n = 0;
    invdiag.clear();
    value.clear();
    rowindex.clear();
    colstart.clear();
    adiag.clear();
  }

  void resize(unsigned int n_) {
    n = n_;
    invdiag.resize(n);
    colstart.resize(n + 1);
    adiag.resize(n);
  }

  void write_matlab(std::ostream &output, const char *variable_name) {
    output << variable_name << "=sparse([";
    for (unsigned int i = 0; i < n; ++i) {
      output << " " << i + 1;
      for (unsigned int j = colstart[i]; j < colstart[i + 1]; ++j) {
        output << " " << rowindex[j] + 1;
      }
    }
    output << "],...\n  [";
    for (unsigned int i = 0; i < n; ++i) {
      output << " " << i + 1;
      for (unsigned int j = colstart[i]; j < colstart[i + 1]; ++j) {
        output << " " << i + 1;
      }
    }
    output << "],...\n  [";
    for (unsigned int i = 0; i < n; ++i) {
      output << " " << (invdiag[i] != 0 ? 1 / invdiag[i] : 0);
      for (unsigned int j = colstart[i]; j < colstart[i + 1]; ++j) {
        output << " " << value[j];
      }
    }
    output << "], " << n << ", " << n << ");" << std::endl;
  }
};

//============================================================================
// Incomplete Cholesky factorization, level zero, with option for modified
// version. Set modification_parameter between zero (regular incomplete
// Cholesky) and one (fully modified version), with values close to one usually
// giving the best results. The min_diagonal_ratio parameter is used to detect
// and correct problems in factorization: if a pivot is this much less than the
// diagonal entry from the original matrix, the original matrix entry is used
// instead.

template <class T>
void factor_modified_incomplete_cholesky0(const SparseMatrix<T> &matrix, SparseColumnLowerFactor<T> &factor, T modification_parameter = 0.97,
                                          T min_diagonal_ratio = 0.25) {
  // first copy lower triangle of matrix into factor (Note: assuming A is
  // symmetric of course!)
  factor.resize(matrix.n);
  zero(factor.invdiag);  // important: eliminate old values from previous solves!
  factor.value.resize(0);
  factor.rowindex.resize(0);
  zero(factor.adiag);
  for (unsigned int i = 0; i < matrix.n; ++i) {
    factor.colstart[i] = (unsigned int)factor.rowindex.size();
    for (unsigned int j = 0; j < matrix.index[i].size(); ++j) {
      if (matrix.index[i][j] > i) {
        factor.rowindex.push_back(matrix.index[i][j]);
        factor.value.push_back(matrix.value[i][j]);
      } else if (matrix.index[i][j] == i) {
        factor.invdiag[i] = factor.adiag[i] = matrix.value[i][j];
      }
    }
  }
  factor.colstart[matrix.n] = (unsigned int)factor.rowindex.size();
  // now do the incomplete factorization (figure out numerical values)

  // MATLAB code:
  // L=tril(A);
  // for k=1:size(L,2)
  //   L(k,k)=sqrt(L(k,k));
  //   L(k+1:end,k)=L(k+1:end,k)/L(k,k);
  //   for j=find(L(:,k))'
  //     if j>k
  //       fullupdate=L(:,k)*L(j,k);
  //       incompleteupdate=fullupdate.*(A(:,j)~=0);
  //       missing=sum(fullupdate-incompleteupdate);
  //       L(j:end,j)=L(j:end,j)-incompleteupdate(j:end);
  //       L(j,j)=L(j,j)-omega*missing;
  //     end
  //   end
  // end

  for (unsigned int k = 0; k < matrix.n; ++k) {
    if (factor.adiag[k] == 0) continue;  // null row/column
    // figure out the final L(k,k) entry
    if (factor.invdiag[k] < min_diagonal_ratio * factor.adiag[k])
      factor.invdiag[k] = 1 / sqrt(factor.adiag[k]);  // drop to Gauss-Seidel here if the pivot
                                                      // looks dangerously small
    else
      factor.invdiag[k] = 1 / sqrt(factor.invdiag[k]);
    // finalize the k'th column L(:,k)
    for (unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p) {
      factor.value[p] *= factor.invdiag[k];
    }
    // incompletely eliminate L(:,k) from future columns, modifying diagonals
    for (unsigned int p = factor.colstart[k]; p < factor.colstart[k + 1]; ++p) {
      unsigned int j = factor.rowindex[p];  // work on column j
      T multiplier = factor.value[p];
      T missing = 0;
      unsigned int a = factor.colstart[k];
      // first look for contributions to missing from dropped entries above the
      // diagonal in column j
      unsigned int b = 0;
      while (a < factor.colstart[k + 1] && factor.rowindex[a] < j) {
        // look for factor.rowindex[a] in matrix.index[j] starting at b
        while (b < matrix.index[j].size()) {
          if (matrix.index[j][b] < factor.rowindex[a])
            ++b;
          else if (matrix.index[j][b] == factor.rowindex[a])
            break;
          else {
            missing += factor.value[a];
            break;
          }
        }
        ++a;
      }
      // adjust the diagonal j,j entry
      if (a < factor.colstart[k + 1] && factor.rowindex[a] == j) {
        factor.invdiag[j] -= multiplier * factor.value[a];
      }
      ++a;
      // and now eliminate from the nonzero entries below the diagonal in column
      // j (or add to missing if we can't)
      b = factor.colstart[j];
      while (a < factor.colstart[k + 1] && b < factor.colstart[j + 1]) {
        if (factor.rowindex[b] < factor.rowindex[a])
          ++b;
        else if (factor.rowindex[b] == factor.rowindex[a]) {
          factor.value[b] -= multiplier * factor.value[a];
          ++a;
          ++b;
        } else {
          missing += factor.value[a];
          ++a;
        }
      }
      // and if there's anything left to do, add it to missing
      while (a < factor.colstart[k + 1]) {
        missing += factor.value[a];
        ++a;
      }
      // and do the final diagonal adjustment from the missing entries
      factor.invdiag[j] -= modification_parameter * multiplier * missing;
    }
  }
}

//============================================================================
// Solution routines with lower triangular matrix.

// solve L*result=rhs
template <class T>
void solve_lower(const SparseColumnLowerFactor<T> &factor, const std::vector<T> &rhs, std::vector<T> &result) {
  assert(factor.n == rhs.size());
  assert(factor.n == result.size());
  result = rhs;
  for (unsigned int i = 0; i < factor.n; ++i) {
    result[i] *= factor.invdiag[i];
    for (unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j) {
      result[factor.rowindex[j]] -= factor.value[j] * result[i];
    }
  }
}

// solve L^T*result=rhs
template <class T>
void solve_lower_transpose_in_place(const SparseColumnLowerFactor<T> &factor, std::vector<T> &x) {
  assert(factor.n == x.size());
  assert(factor.n > 0);
  unsigned int i = factor.n;
  do {
    --i;
    for (unsigned int j = factor.colstart[i]; j < factor.colstart[i + 1]; ++j) {
      x[i] -= factor.value[j] * x[factor.rowindex[j]];
    }
    x[i] *= factor.invdiag[i];
  } while (i != 0);
}

//============================================================================
// Fixed-width triangular rows for row-major five-point grid factors.
template <class T>
struct SparseGridLowerFactor {
  struct Neighbours {
    unsigned int first, second;
    T first_value, second_value;
  };
  std::vector<Neighbours> lower, upper;

  bool construct_from_factor(const SparseColumnLowerFactor<T> &factor) {
    for (unsigned int i = 0; i < factor.n; ++i) {
      if (factor.colstart[i + 1] - factor.colstart[i] > 2) return false;
    }
    lower.assign(factor.n, Neighbours{factor.n, factor.n, T(0), T(0)});
    upper.assign(factor.n, Neighbours{factor.n, factor.n, T(0), T(0)});
    for (unsigned int i = 0; i < factor.n; ++i) {
      for (unsigned int p = factor.colstart[i]; p < factor.colstart[i + 1]; ++p) {
        const unsigned int j = factor.rowindex[p];
        Neighbours &outgoing = upper[i];
        if (outgoing.first == factor.n) {
          outgoing.first = j;
          outgoing.first_value = factor.value[p];
        } else {
          outgoing.second = j;
          outgoing.second_value = factor.value[p];
        }
        // Preserve the source-column order of the generic forward solve.
        Neighbours &incoming = lower[j];
        if (incoming.first == factor.n) {
          incoming.first = i;
          incoming.first_value = factor.value[p];
        } else if (incoming.second == factor.n) {
          incoming.second = i;
          incoming.second_value = factor.value[p];
        } else {
          return false;
        }
      }
    }
    return true;
  }

  void apply(const std::vector<T> &invdiag, const std::vector<T> &rhs, std::vector<T> &result) const {
    const unsigned int n = static_cast<unsigned int>(invdiag.size());
    assert(lower.size() == n && upper.size() == n);
    assert(rhs.size() == n && result.size() == n);
    for (unsigned int i = 0; i < n; ++i) {
      const Neighbours &row = lower[i];
      T value = rhs[i];
      if (row.first != n) value -= row.first_value * result[row.first];
      if (row.second != n) value -= row.second_value * result[row.second];
      result[i] = value * invdiag[i];
    }
    for (unsigned int i = n; i > 0;) {
      --i;
      const Neighbours &row = upper[i];
      T value = result[i];
      if (row.first != n) value -= row.first_value * result[row.first];
      if (row.second != n) value -= row.second_value * result[row.second];
      result[i] = value * invdiag[i];
    }
  }
};

//============================================================================
// Encapsulates the Conjugate Gradient algorithm with incomplete Cholesky
// factorization preconditioner.

template <class T>
struct PCGSolver {
  PCGSolver(void) { set_solver_parameters(1e-4, 100, 0.97, 0.25); }

  // Sub-phase timings (milliseconds) of the most recent solve() call. Only
  // meaningful for profiling; always populated so callers can read them.
  double time_form_preconditioner = 0.0;  // incomplete-Cholesky factorization
  double time_construct_matrix = 0.0;     // dynamic -> fixed-width or CSR conversion
  double time_iterations = 0.0;           // the CG iteration loop
  double time_initial_residual = 0.0;

  // Iteration-loop breakdown, populated only when profile_iterations is enabled.
  double time_matrix_multiply = 0.0;
  double time_apply_preconditioner = 0.0;
  double time_vector_operations = 0.0;

  void set_use_grid_preconditioner(bool enabled) { use_grid_preconditioner = enabled; }
  void set_use_grid_matrix(bool enabled) { use_grid_matrix = enabled; }

  void set_solver_parameters(T tolerance_factor_, int max_iterations_, T modified_incomplete_cholesky_parameter_ = 0.97, T min_diagonal_ratio_ = 0.25) {
    tolerance_factor = tolerance_factor_;
    if (tolerance_factor < 1e-30) tolerance_factor = 1e-30;
    max_iterations = max_iterations_;
    modified_incomplete_cholesky_parameter = modified_incomplete_cholesky_parameter_;
    min_diagonal_ratio = min_diagonal_ratio_;
  }

  bool solve(const SparseMatrix<T> &matrix, const std::vector<T> &rhs, std::vector<T> &result, T &residual_out, int &iterations_out,
             bool profile_iterations = false, bool warm_start = false) {
    using SolverClock = std::chrono::steady_clock;
    auto ms_since = [](std::chrono::time_point<SolverClock> t) {
      return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(SolverClock::now() - t).count()) * 1e-6;
    };
    time_form_preconditioner = time_construct_matrix = time_iterations = time_initial_residual = 0.0;
    time_matrix_multiply = time_apply_preconditioner = time_vector_operations = 0.0;

    unsigned int n = matrix.n;
    if (m.size() != n) {
      m.resize(n);
      s.resize(n);
      z.resize(n);
      r.resize(n);
    }
    auto t_initial_begin = SolverClock::now();
    if (!warm_start || result.size() != n) {
      result.assign(n, T(0));
    } else {
      for (unsigned int i = 0; i < n; ++i) {
        if (matrix.index[i].empty()) result[i] = 0;
      }
    }
    r = rhs;
    const T rhs_norm = BLAS::abs_max(rhs);
    residual_out = rhs_norm;
    if (rhs_norm == 0) {
      zero(result);
      iterations_out = 0;
      time_initial_residual = ms_since(t_initial_begin);
      return true;
    }
    const T tol = tolerance_factor * rhs_norm;
    time_initial_residual = ms_since(t_initial_begin);

    auto t_construct_begin = SolverClock::now();
    grid_matrix_ready = use_grid_matrix && grid_matrix.construct_from_matrix(matrix);
    if (!grid_matrix_ready) fixed_matrix.construct_from_matrix(matrix);
    time_construct_matrix = ms_since(t_construct_begin);
    if (warm_start) {
      auto t_residual_begin = SolverClock::now();
      if (grid_matrix_ready)
        grid_matrix.multiply_and_subtract(result, r);
      else
        multiply_and_subtract(fixed_matrix, result, r);
      residual_out = BLAS::abs_max(r);
      // Discard an unhelpful guess, for example after a large topology change.
      if (residual_out > rhs_norm) {
        zero(result);
        r = rhs;
        residual_out = rhs_norm;
      }
      time_initial_residual += ms_since(t_residual_begin);
    }
    if (residual_out <= tol) {
      iterations_out = 0;
      return true;
    }

    auto t_factor_begin = SolverClock::now();
    form_preconditioner(matrix);
    time_form_preconditioner = ms_since(t_factor_begin);
    apply_preconditioner(r, z);
    scalar rho = BLAS::dot(z, r);
    if (rho == 0 || rho != rho) {
      iterations_out = 0;
      return false;
    }

    s = z;
    auto t_iter_begin = SolverClock::now();
    auto t_kernel_begin = t_iter_begin;
    auto record_kernel = [&](double &elapsed_ms) {
      if (!profile_iterations) return;
      auto now = SolverClock::now();
      elapsed_ms += static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(now - t_kernel_begin).count()) * 1e-6;
      t_kernel_begin = now;
    };
    int iteration;
    for (iteration = 0; iteration < max_iterations; ++iteration) {
      if (grid_matrix_ready)
        grid_matrix.multiply(s, z);
      else
        multiply(fixed_matrix, s, z);
      record_kernel(time_matrix_multiply);
      scalar alpha = rho / BLAS::dot(s, z);
      residual_out = BLAS::update_solution_and_residual(alpha, s, z, result, r);
      record_kernel(time_vector_operations);
      if (residual_out <= tol) {
        iterations_out = iteration + 1;
        time_iterations = ms_since(t_iter_begin);
        return true;
      }
      apply_preconditioner(r, z);
      record_kernel(time_apply_preconditioner);
      scalar rho_new = BLAS::dot(z, r);
      scalar beta = rho_new / rho;
      BLAS::add_scaled(beta, s, z);
      s.swap(z);  // s=beta*s+z
      rho = rho_new;
      record_kernel(time_vector_operations);
    }
    iterations_out = iteration;
    time_iterations = ms_since(t_iter_begin);
    return false;
  }

 protected:
  // internal structures
  SparseColumnLowerFactor<T> ic_factor;  // modified incomplete cholesky factor
  SparseGridLowerFactor<T> grid_factor;
  bool use_grid_preconditioner = false;
  bool grid_preconditioner_ready = false;
  std::vector<T> m, z, s, r;          // temporary vectors for PCG
  FixedSparseMatrix<T> fixed_matrix;  // used within loop
  SparseGridMatrix<T> grid_matrix;
  bool use_grid_matrix = false;
  bool grid_matrix_ready = false;

  // parameters
  T tolerance_factor;
  int max_iterations;
  T modified_incomplete_cholesky_parameter;
  T min_diagonal_ratio;

  void form_preconditioner(const SparseMatrix<T> &matrix) {
    factor_modified_incomplete_cholesky0(matrix, ic_factor);
    grid_preconditioner_ready = use_grid_preconditioner && grid_factor.construct_from_factor(ic_factor);
  }

  void apply_preconditioner(const std::vector<T> &x, std::vector<T> &result) {
    if (grid_preconditioner_ready) {
      grid_factor.apply(ic_factor.invdiag, x, result);
      return;
    }
    solve_lower(ic_factor, x, result);
    solve_lower_transpose_in_place(ic_factor, result);
  }
};
}  // namespace robertbridson

#endif  // APIC2D_PCG_SOLVER_H_
