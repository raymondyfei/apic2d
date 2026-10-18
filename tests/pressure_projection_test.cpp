#include <iostream>
#include <sstream>
#include <stdexcept>

#include "fluidsim.h"

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

struct PressureProjectionTest {
  static constexpr float dt = 0.05f;

  static void initialize(FluidSim& sim) {
    sim.initialize(Vector2s::Zero(), 10.0f, 10, 9, 1.0f, false, false, false, false, true);
    sim.liquid_phi_.assign(0.7f);
    sim.pressure_.assign(sim.ni_ * sim.nj_, 0.0f);
    for (int j = 0; j < sim.u_.nj; ++j) {
      for (int i = 0; i < sim.u_.ni; ++i) {
        sim.u_(i, j) = 0.1f * i - 0.03f * j;
        sim.u_weights_(i, j) = (i + j) % 7 == 0 ? 0.0f : 0.6f;
      }
    }
    for (int j = 0; j < sim.v_.nj; ++j) {
      for (int i = 0; i < sim.v_.ni; ++i) {
        sim.v_(i, j) = -0.05f * i + 0.13f * j;
        sim.v_weights_(i, j) = (i + j) % 5 == 0 ? 0.0f : 0.8f;
      }
    }
    for (int j = 2; j <= 5; ++j) {
      for (int i = 2; i <= 6; ++i) sim.liquid_phi_(i, j) = -0.4f;
    }
    sim.liquid_phi_(2, 2) = -1e-5f;
    sim.liquid_phi_(8, 6) = -0.3f;
    sim.liquid_phi_(1, 3) = -0.4f;
    sim.liquid_phi_(0, 3) = -0.2f;
    sim.u_weights_(3, 3) = sim.u_weights_(4, 3) = 0;
    sim.v_weights_(3, 3) = sim.v_weights_(3, 4) = 0;
  }

  static bool active(const FluidSim& sim, int i, int j) {
    return i > 0 && i < sim.ni_ - 1 && j > 0 && j < sim.nj_ - 1 && sim.liquid_phi_(i, j) < 0 &&
           (sim.u_weights_(i, j) > 0 || sim.u_weights_(i + 1, j) > 0 || sim.v_weights_(i, j) > 0 || sim.v_weights_(i, j + 1) > 0);
  }

  static void assemble_full_grid(const FluidSim& sim, robertbridson::SparseMatrix<float>& matrix, std::vector<float>& rhs, float timestep = dt) {
    matrix.resize(sim.ni_ * sim.nj_);
    matrix.zero();
    rhs.assign(matrix.n, 0.0f);
    const int di[] = {1, -1, 0, 0};
    const int dj[] = {0, 0, 1, -1};
    for (int j = 1; j < sim.nj_ - 1; ++j) {
      for (int i = 1; i < sim.ni_ - 1; ++i) {
        if (!active(sim, i, j)) continue;
        const int row = i + sim.ni_ * j;
        const float phi = sim.liquid_phi_(i, j);
        for (int face = 0; face < 4; ++face) {
          const int fi = i + (di[face] > 0 ? 1 : 0);
          const int fj = j + (dj[face] > 0 ? 1 : 0);
          const float weight = di[face] != 0 ? sim.u_weights_(fi, fj) : sim.v_weights_(fi, fj);
          const float velocity = di[face] != 0 ? sim.u_(fi, fj) : sim.v_(fi, fj);
          const float term = weight * timestep / sqr(sim.dx_);
          const float neighbour_phi = sim.liquid_phi_(i + di[face], j + dj[face]);
          if (neighbour_phi < 0) {
            matrix.add_to_element(row, row, term);
            matrix.add_to_element(row, row + di[face] + sim.ni_ * dj[face], -term);
          } else {
            matrix.add_to_element(row, row, term / std::max(0.01f, phi / (phi - neighbour_phi)));
          }
          rhs[row] -= (di[face] + dj[face]) * weight * velocity / sim.dx_;
        }
      }
    }
  }

  static float theta(float a, float b) {
    if (a < 0 && b < 0) return 1.0f;
    return std::max(0.01f, a < 0 ? a / (a - b) : b / (b - a));
  }

  static void verify_projection(FluidSim& sim, float timestep = dt) {
    robertbridson::SparseMatrix<float> full;
    std::vector<float> rhs;
    assemble_full_grid(sim, full, rhs, timestep);
    std::vector<float> expected_pressure = sim.pressure_;
    robertbridson::PCGSolver<float> reference;
    float residual;
    int iterations;
    require(reference.solve(full, rhs, expected_pressure, residual, iterations, false, sim.pressure_warm_start_), "reference projection failed");
    Array2s expected_u = sim.u_;
    Array2s expected_v = sim.v_;
    for (int j = 0; j < expected_u.nj; ++j) {
      for (int i = 1; i < expected_u.ni - 1; ++i) {
        if (sim.u_weights_(i, j) <= 0) {
          expected_u(i, j) = 0;
        } else if (sim.liquid_phi_(i, j) < 0 || sim.liquid_phi_(i - 1, j) < 0) {
          const int grid = i + sim.ni_ * j;
          expected_u(i, j) -=
              timestep * (expected_pressure[grid] - expected_pressure[grid - 1]) / sim.dx_ / theta(sim.liquid_phi_(i, j), sim.liquid_phi_(i - 1, j));
        }
      }
    }
    for (int j = 1; j < expected_v.nj - 1; ++j) {
      for (int i = 0; i < expected_v.ni; ++i) {
        if (sim.v_weights_(i, j) <= 0) {
          expected_v(i, j) = 0;
        } else if (sim.liquid_phi_(i, j) < 0 || sim.liquid_phi_(i, j - 1) < 0) {
          const int grid = i + sim.ni_ * j;
          expected_v(i, j) -=
              timestep * (expected_pressure[grid] - expected_pressure[grid - sim.ni_]) / sim.dx_ / theta(sim.liquid_phi_(i, j), sim.liquid_phi_(i, j - 1));
        }
      }
    }
    std::ostringstream output;
    auto* old_output = std::cout.rdbuf(output.rdbuf());
    sim.solve_pressure(timestep);
    std::cout.rdbuf(old_output);
    require(output.str().empty(), "compact projection reported a failure");
    int expected_active = 0;
    for (int grid = 0; grid < sim.ni_ * sim.nj_; ++grid) {
      const bool is_active = active(sim, grid % sim.ni_, grid / sim.ni_);
      if (is_active) {
        require(sim.pressure_index_[grid] == expected_active, "compact mapping is not row-major");
        require(sim.pressure_cells_[expected_active] == grid, "compact mapping does not round-trip");
        require(std::fabs(sim.rhs_[expected_active] - rhs[grid]) < 1e-6f, "compact RHS changed");
        ++expected_active;
      } else {
        require(sim.pressure_index_[grid] == -1 && sim.pressure_[grid] == 0, "inactive pressure was not cleared");
      }
      require(std::fabs(sim.pressure_[grid] - expected_pressure[grid]) <= 1e-3f * std::max(1.0f, std::fabs(expected_pressure[grid])),
              "compact pressure differs from full-grid projection");
    }
    require(sim.matrix_.n == expected_active && sim.rhs_.size() == expected_active && sim.active_pressure_.size() == expected_active,
            "solver vectors were not compacted");
    for (int row = 0; row < expected_active; ++row) {
      for (int column = 0; column < expected_active; ++column) {
        const float expected = full(sim.pressure_cells_[row], sim.pressure_cells_[column]);
        require(std::fabs(sim.matrix_(row, column) - expected) < 1e-6f, "compact matrix changed the discretization");
      }
    }
    for (int j = 0; j < expected_u.nj; ++j) {
      for (int i = 0; i < expected_u.ni; ++i) {
        require(std::fabs(sim.u_(i, j) - expected_u(i, j)) <= 1e-3f * std::max(1.0f, std::fabs(expected_u(i, j))), "u projection changed");
      }
    }
    for (int j = 0; j < expected_v.nj; ++j) {
      for (int i = 0; i < expected_v.ni; ++i) {
        require(std::fabs(sim.v_(i, j) - expected_v(i, j)) <= 1e-3f * std::max(1.0f, std::fabs(expected_v(i, j))), "v projection changed");
      }
    }
  }

  static void test_projection_and_empty_sets() {
    for (bool warm : {false, true}) {
      FluidSim sim;
      initialize(sim);
      sim.set_pressure_warm_start(warm);
      verify_projection(sim);
      sim.liquid_phi_.assign(0.7f);
      verify_projection(sim);
      require(sim.matrix_.n == 0, "empty fluid produced unknowns");
      sim.liquid_phi_.assign(-0.4f);
      sim.u_weights_.set_zero();
      sim.v_weights_.set_zero();
      verify_projection(sim);
      require(sim.matrix_.n == 0, "fully blocked fluid produced unknowns");
      sim.liquid_phi_.assign(0.7f);
      for (int j = 1; j < sim.nj_ - 1; ++j) {
        for (int i = 1; i < sim.ni_ - 1; ++i) sim.liquid_phi_(i, j) = -0.4f;
      }
      sim.u_weights_.assign(1.0f);
      sim.v_weights_.assign(1.0f);
      verify_projection(sim);
      require(sim.matrix_.n == (sim.ni_ - 2) * (sim.nj_ - 2), "dense fluid was incorrectly compacted");
    }
  }

  static void test_warm_start_remapping() {
    FluidSim sim;
    initialize(sim);
    verify_projection(sim);
    sim.liquid_phi_(2, 2) = 0.7f;
    sim.liquid_phi_(1, 2) = -0.4f;
    robertbridson::SparseMatrix<float> full;
    std::vector<float> rhs;
    assemble_full_grid(sim, full, rhs);
    std::vector<int> cells;
    for (int grid = 0; grid < sim.ni_ * sim.nj_; ++grid) {
      if (active(sim, grid % sim.ni_, grid / sim.ni_)) cells.push_back(grid);
    }
    const int n = static_cast<int>(cells.size());
    Eigen::MatrixXf dense(n, n);
    Eigen::VectorXf dense_rhs(n);
    for (int row = 0; row < n; ++row) {
      dense_rhs[row] = rhs[cells[row]];
      for (int column = 0; column < n; ++column) dense(row, column) = full(cells[row], cells[column]);
    }
    const Eigen::VectorXf exact = dense.ldlt().solve(dense_rhs);
    sim.pressure_.assign(sim.ni_ * sim.nj_, 100.0f);
    for (int row = 0; row < n; ++row) sim.pressure_[cells[row]] = exact[row];
    const long long previous_iterations = sim.timing_pcg_iterations_;
    sim.solver_.set_solver_parameters(1e-4f, 0);
    std::ostringstream output;
    auto* old_output = std::cout.rdbuf(output.rdbuf());
    sim.solve_pressure(dt);
    std::cout.rdbuf(old_output);
    require(output.str().empty() && sim.timing_pcg_iterations_ == previous_iterations, "warm start used obsolete compact indices");
    for (int row = 0; row < n; ++row) require(sim.pressure_[cells[row]] == exact[row], "remapped warm-start pressure changed");
    require(sim.pressure_[2 + sim.ni_ * 2] == 0.0f, "newly inactive pressure survived scattering");

    sim.timing_step_count_ = 1;
    output.str("");
    old_output = std::cout.rdbuf(output.rdbuf());
    sim.report_timing();
    sim.report_timing();
    std::cout.rdbuf(old_output);
    require(output.str().find("active unknowns:") != std::string::npos && output.str().find("[pressure map]") != std::string::npos &&
                output.str().find("[pressure scatter]") != std::string::npos && output.str().find("[pressure solve total]") != std::string::npos &&
                output.str().find("assembly row reuse:") != std::string::npos,
            "compaction profiling is incomplete");
    const auto summary = output.str().find("[pressure solve total]");
    require(output.str().find("[pressure solve total]", summary + 1) == std::string::npos, "inclusive pressure timing was printed more than once");
    require(sim.timing_pressure_active_cells_ == 0 && sim.timing_pressure_grid_cells_ == 0 && sim.timing_pressure_solve_ == 0 &&
                sim.timing_pressure_reused_rows_ == 0,
            "compaction profiling was not reset");
  }

  static void test_row_structure_reuse() {
    FluidSim sim;
    initialize(sim);
    sim.set_pressure_warm_start(false);
    verify_projection(sim);
    require(sim.timing_pressure_reused_rows_ == 0, "first projection reused uninitialized assembly rows");
    sim.u_(4, 2) += 0.3f;
    sim.v_(2, 4) -= 0.2f;
    sim.u_weights_(4, 2) *= 0.5f;
    sim.liquid_phi_(2, 2) = -0.2f;
    const auto previous_rows = sim.timing_pressure_reused_rows_;
    verify_projection(sim, 0.5f * dt);
    require(sim.timing_pressure_reused_rows_ == previous_rows + sim.matrix_.n, "matching assembly rows were not reused");
    const unsigned int active_count = sim.matrix_.n;
    const auto before_removal = sim.timing_pressure_reused_rows_;
    sim.u_weights_(4, 2) = 0.0f;
    sim.v_(2, 4) -= 0.2f;
    verify_projection(sim);
    require(sim.matrix_.n == active_count && sim.timing_pressure_reused_rows_ < before_removal + active_count,
            "face weight crossing zero did not rebuild affected rows");
    const auto before_restoration = sim.timing_pressure_reused_rows_;
    sim.u_weights_(4, 2) = 0.6f;
    sim.v_(2, 4) -= 0.2f;
    verify_projection(sim);
    require(sim.timing_pressure_reused_rows_ < before_restoration + active_count, "restored connectivity reused stale row indices");
  }

  static void test_disabled_profiling() {
    FluidSim sim;
    initialize(sim);
    sim.print_timing_ = false;
    verify_projection(sim);
    require(sim.timing_accum_.empty() && sim.timing_pressure_active_cells_ == 0 && sim.timing_pressure_grid_cells_ == 0 && sim.timing_pressure_solve_ == 0 &&
                sim.timing_pressure_reused_rows_ == 0,
            "disabled profiling recorded compaction statistics");
  }

  static void extrapolate_reference(Array2s& grid, const Array2s& weights, const Array2s& phi, const Vector2i& offset, Array2c& valid) {
    valid.resize(grid.ni, grid.nj);
    valid.set_zero();
    for (int j = 1; j < grid.nj - 1; ++j) {
      for (int i = 1; i < grid.ni - 1; ++i) valid(i, j) = weights(i, j) > 0 && (phi(i, j) < 0 || phi(i + offset[0], j + offset[1]) < 0);
    }
    for (int layer = 0; layer < 2; ++layer) {
      const Array2s previous = grid;
      const Array2c known = valid;
      for (int j = 1; j < grid.nj - 1; ++j) {
        for (int i = 1; i < grid.ni - 1; ++i) {
          if (known(i, j)) continue;
          float sum = 0;
          int count = 0;
          for (const Vector2i& neighbour : {Vector2i(i + 1, j), Vector2i(i - 1, j), Vector2i(i, j + 1), Vector2i(i, j - 1)}) {
            if (known(neighbour[0], neighbour[1])) {
              sum += previous(neighbour[0], neighbour[1]);
              ++count;
            }
          }
          if (count) {
            grid(i, j) = sum / count;
            valid(i, j) = 1;
          }
        }
      }
    }
  }

  static void test_frontier_extrapolation() {
    for (int width : {2, 3, 13}) {
      FluidSim sim;
      sim.initialize(Vector2s::Zero(), 10, width, width - 1, 1, false, false, false, false, false);
      for (int pattern = 0; pattern < 6; ++pattern) {
        for (int j = 0; j < sim.liquid_phi_.nj; ++j) {
          for (int i = 0; i < sim.liquid_phi_.ni; ++i) {
            const bool fluid = pattern == 1 || pattern == 5 || (pattern == 2 && i == width / 2 && j == width / 2) || (pattern == 3 && (i + 3 * j) % 5 == 0) ||
                               (pattern == 4 && (i < 2 || j < 2));
            sim.liquid_phi_(i, j) = fluid ? -0.5f : 0.5f;
          }
        }
        for (bool horizontal : {true, false}) {
          Array2s& grid = horizontal ? sim.u_ : sim.v_;
          Array2s& weights = horizontal ? sim.u_weights_ : sim.v_weights_;
          const Vector2i offset = horizontal ? Vector2i(-1, 0) : Vector2i(0, -1);
          for (int j = 0; j < grid.nj; ++j) {
            for (int i = 0; i < grid.ni; ++i) {
              grid(i, j) = ((i * 7 + j * 3) % 17 - 8) * 0.125f;
              weights(i, j) = pattern == 5 || (pattern == 3 && (i + j) % 3 == 0) ? 0 : 0.7f;
            }
          }
          Array2s expected = grid;
          Array2c expected_valid;
          extrapolate_reference(expected, weights, sim.liquid_phi_, offset, expected_valid);
          sim.valid_.assign(1);
          sim.extrapolate_velocity(grid, weights, offset);
          require(sim.valid_.ni == expected_valid.ni && sim.valid_.nj == expected_valid.nj, "validity scratch dimensions did not match the component");
          for (int j = 0; j < grid.nj; ++j) {
            for (int i = 0; i < grid.ni; ++i) {
              require(sim.valid_(i, j) == expected_valid(i, j), "frontier propagation leaked across layers, components, or boundaries");
              require(std::fabs(grid(i, j) - expected(i, j)) < 1e-6f, "frontier extrapolation changed Jacobi averaging");
            }
          }
        }
      }
    }
  }

  static void test_large_grid_projection() {
    FluidSim sim;
    sim.initialize(Vector2s::Zero(), 100, 256, 256, 1, false, false, false, false, false);
    sim.liquid_phi_.assign(0.7f);
    sim.u_weights_.assign(0.8f);
    sim.v_weights_.assign(0.6f);
    for (int j = 0; j < sim.u_.nj; ++j) {
      for (int i = 0; i < sim.u_.ni; ++i) sim.u_(i, j) = 0.01f * ((i + j) % 7);
    }
    for (int j = 0; j < sim.v_.nj; ++j) {
      for (int i = 0; i < sim.v_.ni; ++i) sim.v_(i, j) = -0.01f * ((2 * i + j) % 11);
    }
    for (int j = 100; j < 106; ++j) {
      for (int i = 100; i < 106; ++i) sim.liquid_phi_(i, j) = -0.4f;
    }
    verify_projection(sim);
  }
};

int main() {
  PressureProjectionTest::test_projection_and_empty_sets();
  PressureProjectionTest::test_warm_start_remapping();
  PressureProjectionTest::test_row_structure_reuse();
  PressureProjectionTest::test_disabled_profiling();
  PressureProjectionTest::test_frontier_extrapolation();
  PressureProjectionTest::test_large_grid_projection();
  std::cout << "Compact pressure projection and extrapolation checks passed\n";
}
