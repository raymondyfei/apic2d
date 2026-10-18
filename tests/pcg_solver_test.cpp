#include "pcg_solver.h"

#include <iostream>
#include <limits>
#include <stdexcept>

using robertbridson::PCGSolver;
using robertbridson::SparseMatrix;

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

struct PreconditionerTestSolver : PCGSolver<float> {
  using PCGSolver<float>::apply_preconditioner;
  using PCGSolver<float>::form_preconditioner;
  using PCGSolver<float>::grid_preconditioner_ready;
  using PCGSolver<float>::ic_factor;
  using PCGSolver<float>::grid_matrix_ready;
};

template <class T>
void test_grid_matrix_products() {
  robertbridson::SparseGridMatrix<T> stencil;
  for (unsigned int n : {8, 1, 0}) {
    SparseMatrix<T> matrix(n);
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < std::min(n, i % 6); ++j) matrix.set_element(i, (i + j) % n, T(0.125) * (i + 1) - T(0.25) * j);
    }
    for (int pass = 0; pass < 2; ++pass) {
      if (pass && n) matrix.set_element(n - 1, 0, T(-0.75));
      require(stencil.construct_from_matrix(matrix), "five-entry matrix was rejected");
      robertbridson::FixedSparseMatrix<T> csr;
      csr.construct_from_matrix(matrix);
      std::vector<T> x(n), actual(n), expected(n);
      for (unsigned int i = 0; i < n; ++i) x[i] = T(0.5) * i - T(1);
      robertbridson::multiply(csr, x, expected);
      stencil.multiply(x, actual);
      T scale = 1;
      for (T value : expected) scale = std::max(scale, std::fabs(value));
      const T tolerance = 32 * std::numeric_limits<T>::epsilon() * scale;
      for (unsigned int i = 0; i < n; ++i) require(std::fabs(actual[i] - expected[i]) <= tolerance, "fixed-width product differs from CSR");
      actual.assign(n, T(0.25));
      expected = actual;
      robertbridson::multiply_and_subtract(csr, x, expected);
      stencil.multiply_and_subtract(x, actual);
      for (unsigned int i = 0; i < n; ++i) require(std::fabs(actual[i] - expected[i]) <= tolerance, "fixed-width warm residual differs from CSR");
    }
  }
}

void test_grid_matrix_fallback() {
  PreconditionerTestSolver specialized;
  specialized.set_use_grid_matrix(true);
  PCGSolver<float> generic;
  SparseMatrix<float> matrix(6);
  std::vector<float> rhs(6, 1.0f), actual(6), expected(6);
  for (int pass = 0; pass < 4; ++pass) {
    const bool eligible = pass == 1;
    matrix.zero();
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        if (pass == 1 || pass == 3) {
          if (i == j) matrix.set_element(i, j, 2.0f + i);
        } else {
          matrix.set_element(i, j, i == j ? 3.0f : -0.125f);
        }
      }
    }
    specialized.set_use_grid_matrix(pass != 3);
    for (bool warm : {false, true}) {
      rhs[0] += 0.01f;
      float residual, reference_residual;
      int iterations, reference_iterations;
      require(generic.solve(matrix, rhs, expected, reference_residual, reference_iterations, false, warm), "reference matrix solve failed");
      require(specialized.solve(matrix, rhs, actual, residual, iterations, false, warm), "grid matrix or fallback solve failed");
      require(specialized.grid_matrix_ready == eligible, "grid matrix dispatch did not follow eligibility or the opt-in setting");
      for (int i = 0; i < 6; ++i) require(std::fabs(actual[i] - expected[i]) < 3e-4f, "grid matrix fallback changed the solution");
    }
  }
}

void compare_preconditioner_application(PreconditionerTestSolver& solver, const SparseMatrix<float>& matrix, bool eligible) {
  solver.form_preconditioner(matrix);
  require(solver.grid_preconditioner_ready == eligible, "grid preconditioner eligibility is incorrect");
  std::vector<float> rhs(matrix.n), expected(matrix.n), actual(matrix.n);
  for (unsigned int i = 0; i < matrix.n; ++i) rhs[i] = static_cast<float>(static_cast<int>(i % 7) - 3) * 0.25f;
  robertbridson::solve_lower(solver.ic_factor, rhs, expected);
  robertbridson::solve_lower_transpose_in_place(solver.ic_factor, expected);
  solver.apply_preconditioner(rhs, actual);
  const float tolerance = 2e-6f * std::max(1.0f, robertbridson::BLAS::abs_max(expected));
  for (unsigned int i = 0; i < matrix.n; ++i) {
    require(std::fabs(actual[i] - expected[i]) <= tolerance, "specialized MIC application differs from generic solves");
  }
}

void test_grid_preconditioner() {
  PreconditionerTestSolver specialized;
  specialized.set_use_grid_preconditioner(true);
  specialized.set_use_grid_matrix(true);
  for (int width : {1, 2, 7, 12}) {
    std::vector<int> cells(width * width, -1);
    int count = 0;
    for (int grid = 0; grid < width * width; ++grid) {
      if (grid % 5 != 2) cells[grid] = count++;
    }
    SparseMatrix<float> matrix(count + 1);
    for (int grid = 0; grid < width * width; ++grid) {
      const int row = cells[grid];
      if (row < 0) continue;
      matrix.set_element(row, row, 4.5f + 0.1f * (grid % 3));
      for (int offset : {1, width}) {
        if (offset == 1 && grid % width == width - 1) continue;
        const int neighbour = grid + offset;
        if (neighbour >= width * width || cells[neighbour] < 0) continue;
        const float value = -0.2f - 0.1f * (grid % 3);
        matrix.set_element(row, cells[neighbour], value);
        matrix.set_element(cells[neighbour], row, value);
      }
    }
    compare_preconditioner_application(specialized, matrix, true);
    std::vector<float> rhs(matrix.n, 1.0f), generic_result(matrix.n), grid_result(matrix.n);
    rhs.back() = 0.0f;
    PCGSolver<float> generic;
    for (bool warm : {false, true}) {
      rhs[0] += 0.001f;
      float generic_residual, grid_residual;
      int generic_iterations, grid_iterations;
      require(generic.solve(matrix, rhs, generic_result, generic_residual, generic_iterations, false, warm), "generic MIC solve failed");
      require(specialized.solve(matrix, rhs, grid_result, grid_residual, grid_iterations, true, warm), "specialized MIC solve failed");
      require(grid_iterations <= generic_iterations + 1, "specialized application unexpectedly increased iterations");
      for (unsigned int i = 0; i < matrix.n; ++i) {
        require(std::fabs(grid_result[i] - generic_result[i]) <= 3e-4f * std::max(1.0f, std::fabs(generic_result[i])),
                "specialized application changed the PCG solution");
      }
    }
  }

  for (int centre : {0, 3}) {
    SparseMatrix<float> matrix(4);
    for (int i = 0; i < 4; ++i) {
      matrix.set_element(i, i, 4.0f);
      if (i != centre) {
        matrix.set_element(i, centre, -0.2f);
        matrix.set_element(centre, i, -0.2f);
      }
    }
    compare_preconditioner_application(specialized, matrix, false);
  }
  SparseMatrix<float> diagonal(2);
  diagonal.set_element(0, 0, 2.0f);
  diagonal.set_element(1, 1, 3.0f);
  compare_preconditioner_application(specialized, diagonal, true);
  specialized.set_use_grid_preconditioner(false);
  compare_preconditioner_application(specialized, diagonal, false);

  robertbridson::SparseColumnLowerFactor<float> empty;
  robertbridson::SparseGridLowerFactor<float> empty_grid;
  require(empty_grid.construct_from_factor(empty), "empty grid factor construction failed");
  std::vector<float> no_values;
  empty_grid.apply(empty.invdiag, no_values, no_values);
}

void test_fused_update() {
  for (size_t size : {0, 1, 3, 4, 5, 7, 16, 17, 31, 32, 33, 127, 10003}) {
    std::vector<float> direction(size), product(size), solution(size), residual(size);
    for (size_t i = 0; i < size; ++i) {
      direction[i] = static_cast<float>(static_cast<int>(i % 19) - 9) * 0.125f;
      product[i] = static_cast<float>(static_cast<int>(i % 23) - 11) * 0.25f;
      solution[i] = static_cast<float>(i % 13) * 0.5f;
      residual[i] = static_cast<float>(static_cast<int>(i % 17) - 8) * 0.125f;
    }
    std::vector<float> expected_solution = solution, expected_residual = residual;
    robertbridson::BLAS::add_scaled(0.125f, direction, expected_solution);
    robertbridson::BLAS::add_scaled(-0.125f, product, expected_residual);
    float expected_norm = 0.0f;
    for (float value : expected_residual) expected_norm = std::max(expected_norm, std::fabs(value));
    const float norm = robertbridson::BLAS::update_solution_and_residual(0.125f, direction, product, solution, residual);
    require(solution == expected_solution && residual == expected_residual, "fused update differs from separate updates");
    require(norm == expected_norm && robertbridson::BLAS::abs_max(residual) == expected_norm, "fused or standalone norm is incorrect for a vector tail");
  }
}

void test_warm_start(bool grid_matrix) {
  const unsigned int width = 12;
  const unsigned int n = width * width;
  SparseMatrix<float> matrix(n);
  for (unsigned int i = 0; i < n; ++i) {
    matrix.set_element(i, i, 4.5f);
    if (i % width > 0) matrix.set_element(i, i - 1, -1.0f);
    if (i % width + 1 < width) matrix.set_element(i, i + 1, -1.0f);
    if (i >= width) matrix.set_element(i, i - width, -1.0f);
    if (i + width < n) matrix.set_element(i, i + width, -1.0f);
  }

  PCGSolver<float> solver;
  solver.set_use_grid_matrix(grid_matrix);
  std::vector<float> rhs(n, 1.0f), previous(n), cold(n);
  float residual;
  int iterations;
  require(solver.solve(matrix, rhs, previous, residual, iterations), "initial cold solve failed");
  rhs[n / 2] += 0.001f;
  require(solver.solve(matrix, rhs, cold, residual, iterations), "updated cold solve failed");
  const int cold_iterations = iterations;
  require(solver.solve(matrix, rhs, previous, residual, iterations, true, true), "warm solve failed");
  require(iterations < cold_iterations, "nearby-system warm start did not reduce iterations");
  std::vector<float> true_residual = rhs;
  robertbridson::multiply_and_subtract(matrix, previous, true_residual);
  require(robertbridson::BLAS::abs_max(true_residual) <= 1.1e-4f * robertbridson::BLAS::abs_max(rhs), "warm solve violated the RHS-relative tolerance");
  require(solver.time_initial_residual > 0.0, "initial residual calculation was not timed");

  previous.assign(n, 1e6f);
  require(solver.solve(matrix, rhs, previous, residual, iterations, false, true), "bad-guess fallback failed");
  require(iterations == cold_iterations && previous == cold, "bad guess did not fall back to the cold solve");
}

void test_initial_guess_and_residual_norm() {
  PCGSolver<float> solver;
  SparseMatrix<float> matrix(2);
  matrix.set_element(0, 0, 1.0f);
  matrix.set_element(1, 1, 1.0f);
  std::vector<float> rhs{1.0f, 1.0f}, result{1.0f, 1.0f};
  float residual;
  int iterations;
  require(solver.solve(matrix, rhs, result, residual, iterations, true, true), "exact-guess solve failed");
  require(iterations == 0 && residual == 0 && solver.time_form_preconditioner == 0.0 && solver.time_iterations == 0.0, "exact guess did not skip PCG");
  require(solver.solve(matrix, rhs, result, residual, iterations), "default cold solve failed");
  require(iterations == 1 && result == rhs, "cold solve reused the previous pressure");

  solver.set_solver_parameters(1e-4f, 0);
  result = {1.00005f, 1.0f};
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "tolerance was based on the initial residual instead of RHS");
  require(iterations == 0 && result[0] > 1.0f, "acceptable initial guess was discarded");

  solver.set_solver_parameters(1e-4f, 100);
  result = {1.5f, 1.0f};
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "negative-residual solve failed");
  require(iterations == 1 && result == rhs, "negative residual incorrectly passed convergence check");
  require(robertbridson::BLAS::abs_max(std::vector<float>{-10.0f, 1.0f, 0.0f}) == 10.0f, "mixed-sign infinity norm is incorrect");

  result = {0.0f, 0.0f};
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "zero-guess warm solve failed");
  require(iterations == 1 && result == rhs, "zero guess did not behave like cold start");
  rhs.assign(2, 0.0f);
  require(solver.solve(matrix, rhs, result, residual, iterations, true, true), "zero-RHS solve failed");
  require(iterations == 0 && result == rhs && solver.time_iterations == 0.0 && solver.time_matrix_multiply == 0.0,
          "zero RHS retained stale pressure or timings");
}

void test_inactive_cells_and_resizing() {
  PCGSolver<float> solver;
  SparseMatrix<float> matrix(2);
  matrix.set_element(0, 0, 1.0f);
  std::vector<float> rhs{1.0f, 0.0f}, result{1.0f, 100.0f};
  float residual;
  int iterations;
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "inactive-cell solve failed");
  require(iterations == 0 && result[1] == 0.0f, "inactive pressure was not cleared");
  matrix.set_element(1, 1, 1.0f);
  rhs[1] = 2.0f;
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "newly active-cell solve failed");
  require(result == rhs, "newly active pressure was not solved");

  result.resize(1);
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "resized-guess solve failed");
  require(result == rhs, "mismatched initial guess size was not handled");
  matrix.resize(0);
  rhs.clear();
  require(solver.solve(matrix, rhs, result, residual, iterations, false, true), "empty-system solve failed");
  require(result.empty() && iterations == 0 && residual == 0, "empty-system result is incorrect");
}

int main() {
  test_grid_matrix_products<float>();
  test_grid_matrix_products<double>();
  test_grid_matrix_fallback();
  test_grid_preconditioner();
  test_fused_update();
  test_warm_start(false);
  test_warm_start(true);
  test_initial_guess_and_residual_norm();
  test_inactive_cells_and_resizing();
  std::cout << "PCG solver regression checks passed\n";
}
