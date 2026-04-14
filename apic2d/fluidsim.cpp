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

#ifdef WIN32
#define NOMINMAX
#endif

#include "fluidsim.h"

#include "array2_utils.h"
#include "kernel.h"
#include "sorter.h"
#include "sparse_matrix.h"

#ifdef __APPLE__
#include <GLUT/glut.h>  // why does Apple have to put glut.h here...
#else
#include <GL/glut.h>  // ...when everyone else puts it here?
#endif

#include "openglutils.h"

// Change here to try different integration scheme, options:
// IT_PIC: original particle-in-cell (PIC)
// IT_FLIP: original fluid-implicit-particle (FLIP)
// IT_RPIC: rotational particle-in-cell (RPIC)
// IT_APIC: affine particle-in-cell (APIC)
// IT_AFLIP: affine fluid-implicit-particle (AFLIP)
// IT_ASFLIP: affine separable fluid-implicit-particle (ASFLIP)
const FluidSim::INTEGRATOR_TYPE integration_scheme = FluidSim::IT_APIC;

// Change here to try different order for velocity evaluation from grid,
// options:
// VO_EULER: first order evaluation
// VO_RA2: Ralston's second order
//         evaluation
// VO_RK3: Runge Kutta's 3rd-order method
// VO_RK4: Runge Kutta's
//         4rd-order method
const FluidSim::VELOCITY_ORDER velocity_order = FluidSim::VO_EULER;

// Change here to try different order for interpolation
// options:
// IO_LINEAR: linear interpolation
// IO_QUADRATIC: quadratic interpolation
const FluidSim::INTERPOLATION_ORDER interpolation_order = FluidSim::IO_LINEAR;

const scalar lagrangian_ratio = 0.97f;
const int particle_correction_step = 1;

// whether collision is resolved with another pass on each particle (instead of
// grid)
const bool use_particle_boundary_collisions = false;

scalar fraction_inside(scalar phi_left, scalar phi_right);
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid_, Array2c old_valid_,
                 const Vector2i& offset, int num_layers);

FluidSim::~FluidSim() { delete m_sorter_; }

void FluidSim::initialize(const Vector2s& origin, scalar width, int ni, int nj, scalar rho, bool draw_grid, bool draw_particles, bool draw_velocities,
                          bool draw_boundaries, bool print_timing) {
  rho_ = rho;
  draw_grid_ = draw_grid;
  draw_particles_ = draw_particles;
  draw_velocities_ = draw_velocities;
  draw_boundaries_ = draw_boundaries;
  print_timing_ = print_timing;
  origin_ = origin;
  ni_ = ni;
  nj_ = nj;
  dx_ = width / static_cast<scalar>(ni_);
  u_.resize(ni_ + 1, nj_);
  temp_u_.resize(ni_ + 1, nj_);
  u_weights_.resize(ni_ + 1, nj_);
  u_valid_.resize(ni_ + 1, nj_);
  saved_u_.resize(ni_ + 1, nj_);
  v_.resize(ni_, nj_ + 1);
  temp_v_.resize(ni_, nj_ + 1);
  v_weights_.resize(ni_, nj_ + 1);
  v_valid_.resize(ni_, nj_ + 1);
  saved_v_.resize(ni_, nj_ + 1);
  u_.set_zero();
  v_.set_zero();
  nodal_solid_phi_.resize(ni_ + 1, nj_ + 1);
  valid_.resize(ni_ + 1, nj_ + 1);
  old_valid_.resize(ni_ + 1, nj_ + 1);
  liquid_phi_.resize(ni_, nj_);
  m_sorter_ = new sorter(ni_, nj_);
}

/*!
  \brief  Initialize the grid-based signed distance field that dictates the position of
          the solid boundary
*/
void FluidSim::update_boundary() {
  parallel_for(0, nj_ + 1, [this](int j) {
    for (int i = 0; i < ni_ + 1; ++i) {
      Vector2s pos(i * dx_, j * dx_);
      nodal_solid_phi_(i, j) = solid_distance(pos + origin_);
    }
  });
}

void FluidSim::paint_velocity(const Vector2s& brush_center, const scalar brush_radius, const Vector2s& vel) {
  velocity_brush_data_.emplace_back(brush_center, vel, brush_radius);
}

void FluidSim::relaxation(scalar dt) {
  int np = static_cast<int>(particles_.size());

  const scalar re = dx_ / sqrt(2.0) * 1.1;

  int offset = rand() % particle_correction_step;
  // Compute Pseudo Moved Point
  parallel_for(0, np, [&](int n) {
    if (n % particle_correction_step != offset) return;

    Particle& p = particles_[n];

    Vector2s spring = Vector2s::Zero();

    int ix = std::max(0, std::min(static_cast<int>((p.x_(0) - origin_(0)) / dx_), ni_));
    int iy = std::max(0, std::min(static_cast<int>((p.x_(1) - origin_(1)) / dx_), nj_));

    m_sorter_->getNeigboringParticles_cell(ix, iy, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
      for (const Particle* np : neighbors) {
        if (&p != np) {
          scalar dist = (p.x_ - np->x_).norm();
          scalar w = 50.0 * kernel::smooth_kernel(dist * dist, re);
          if (dist > 0.01 * re) {
            spring += w * (p.x_ - np->x_) / dist * re;
          } else {
            spring(0) += 0.01 * re / dt * (rand() & 0xFF) / 255.0;
            spring(1) += 0.01 * re / dt * (rand() & 0xFF) / 255.0;
          }
        }
      }
    });

    p.buf0_ = p.x_ + dt * spring;

    Vector2s pp = (p.buf0_ - origin_) / dx_;
    scalar phi_value = interpolate_value(pp, nodal_solid_phi_);

    if (phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi_);
      normal.normalize();
      p.buf0_ -= phi_value * normal;
    }
  });

  // Update
  parallel_for(0, np, [&](int n) {
    if (n % particle_correction_step != offset) return;
    Particle& p = particles_[n];

    p.x_ = p.buf0_;
  });
}

/*!
  \brief  The main fluid simulation step
*/
void FluidSim::advance(scalar dt) {
  // Passively advect particles_
  tick();
  m_sorter_->sort(particles_, origin_, dx_);
  tock("sort");

  tick();
  map_p2g();
  tock("p2g");
  save_velocity();

  tick();
  add_force(dt);
  tock("add force");

  tick();
  compute_liquid_distance();
  tock("compute phi");

  // Compute finite-volume type_ face area weight for each velocity sample.
  tick();
  compute_weights();
  tock("compute weights");

  // Set up and solve the variational pressure solve.
  tick();
  solve_pressure(dt);
  tock("solve pressure");

  // Pressure projection only produces valid velocities in faces with non-zero
  // associated face area. Because the advection step may interpolate from these
  // invalid faces, we must extrapolate velocities from the fluid domain into
  // these zero-area faces.
  tick();
  extrapolate(u_, temp_u_, u_weights_, liquid_phi_, valid_, old_valid_, Vector2i(-1, 0), 2);
  extrapolate(v_, temp_v_, v_weights_, liquid_phi_, valid_, old_valid_, Vector2i(0, -1), 2);
  tock("extrapolate");

  // For extrapolated velocities, replace the normal component with
  // that of the object.
  tick();
  constrain_velocity();
  tock("constrain velocity");

  tick();
  relaxation(dt);
  tock("relaxation");

  tick();
  switch (integration_scheme) {
    case IT_PIC:
      // The original PIC, which is a specific case of a more general (Affine)
      // FLIP scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 0.0, 0.0);
      break;

    case IT_FLIP:
      // FLIP scheme
      map_g2p_flip_general(dt, lagrangian_ratio, 0.0, 0.0, 0.0);
      break;

    case IT_RPIC:
      // RPIC is a rotation-augmented PIC scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 0.0, 1.0);
      break;

    case IT_APIC:
      // APIC is an affine-augmented PIC scheme
      map_g2p_flip_general(dt, 0.0, 0.0, 1.0, 1.0);
      break;

    case IT_AFLIP:
      // AFLIP is an affine-augmented FLIP scheme
      map_g2p_flip_general(dt, lagrangian_ratio, 0.0, 1.0, 1.0);
      break;

    case IT_ASFLIP:
      // ASFLIP is an affine-augmented PIC scheme with easier particle
      // separation
      map_g2p_flip_general(dt, lagrangian_ratio, 1.0, 1.0, 1.0);
      break;

    default:
      std::cerr << "Unknown integrator type_!" << std::endl;
      break;
  }
  tock("g2p");

  if (use_particle_boundary_collisions) {
    tick();
    particle_boundary_collision(dt);
    tock("particle collision");
  }
}

void FluidSim::save_velocity() {
  if (lagrangian_ratio > 0.0) {
    saved_u_ = u_;
    saved_v_ = v_;
  }
}

void FluidSim::add_force(scalar dt) {
  // gravity
  for (int j = 0; j < nj_ + 1; ++j) {
    for (int i = 0; i < ni_; ++i) {
      v_(i, j) += -981.0 * dt;
    }
  }

  // user-specified velocity
  for (const BrushFootprint& p : velocity_brush_data_) {
    scalar radius_in_grid = p.radius_ / dx_;
    Vector2s center_in_grid_u = (p.center_ - origin_) / dx_ - Vector2s(0.0, 0.5);
    Vector2i low_range_u = Vector2i(std::max(0, std::min(u_.ni - 1, static_cast<int>(floor(center_in_grid_u(0) - radius_in_grid)))),
                                    std::max(0, std::min(u_.nj - 1, static_cast<int>(floor(center_in_grid_u(1) - radius_in_grid)))));
    Vector2i high_range_u = Vector2i(std::max(0, std::min(u_.ni - 1, static_cast<int>(ceil(center_in_grid_u(0) + radius_in_grid)))),
                                     std::max(0, std::min(u_.nj - 1, static_cast<int>(ceil(center_in_grid_u(1) + radius_in_grid)))));
    for (int i = low_range_u(0); i <= high_range_u(0); ++i) {
      for (int j = low_range_u(1); j <= high_range_u(1); ++j) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        if ((pos - p.center_).norm() <= p.radius_) {
          u_(i, j) = p.vel_(0);
        }
      }
    }

    Vector2s center_in_grid_v = (p.center_ - origin_) / dx_ - Vector2s(0.5, 0.0);
    Vector2i low_range_v = Vector2i(std::max(0, std::min(v_.ni - 1, static_cast<int>(floor(center_in_grid_v(0) - radius_in_grid)))),
                                    std::max(0, std::min(v_.nj - 1, static_cast<int>(floor(center_in_grid_v(1) - radius_in_grid)))));
    Vector2i high_range_v = Vector2i(std::max(0, std::min(v_.ni - 1, static_cast<int>(ceil(center_in_grid_v(0) + radius_in_grid)))),
                                     std::max(0, std::min(v_.nj - 1, static_cast<int>(ceil(center_in_grid_v(1) + radius_in_grid)))));
    for (int i = low_range_v(0); i <= high_range_v(0); ++i) {
      for (int j = low_range_v(1); j <= high_range_v(1); ++j) {
        Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
        if ((pos - p.center_).norm() <= p.radius_) {
          v_(i, j) = p.vel_(1);
        }
      }
    }
  }

  velocity_brush_data_.clear();
}

scalar FluidSim::compute_cfl() const {
  scalar max_vel_squared = 0.0;
  for (const Particle& p : particles_) {
    max_vel_squared = std::max(max_vel_squared, p.v_.squaredNorm());
  }
  return dx_ / sqrt(max_vel_squared);
}

/*!
  \brief  For extrapolated points, replace the normal component
          of velocity with the object velocity (in this case zero).
*/
void FluidSim::constrain_velocity() {
  // At lower grid resolutions, the normal estimate from the signed
  // distance function is poor, so it doesn't work quite as well.
  // An exact normal would do better.

  // constrain u
  parallel_for(0, std::max(u_.nj, v_.nj), [this](int j) {
    if (j < u_.nj) {
      for (int i = 0; i < u_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
          Vector2s vel = get_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i, j + 0.5), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_u_(i, j) = vel[0] - perp_component * normal[0];
        }
      }
    }

    if (j < v_.nj) {
      for (int i = 0; i < v_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          // apply constraint
          Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
          Vector2s vel = get_velocity(pos);
          Vector2s normal(0, 0);
          interpolate_gradient(normal, Vector2s(i + 0.5, j), nodal_solid_phi_);
          normal.normalize();
          scalar perp_component = vel.dot(normal);
          temp_v_(i, j) = vel[1] - perp_component * normal[1];
        }
      }
    }
  });

  parallel_for(0, std::max(u_.nj, v_.nj), [this](int j) {
    if (j < u_.nj) {
      for (int i = 0; i < u_.ni; ++i) {
        if (u_weights_(i, j) == 0) {
          u_(i, j) = temp_u_(i, j);
        }
      }
    }
    if (j < v_.nj) {
      for (int i = 0; i < v_.ni; ++i) {
        if (v_weights_(i, j) == 0) {
          v_(i, j) = temp_v_(i, j);
        }
      }    
    }
  });
}

void FluidSim::compute_liquid_distance() {
  const scalar min_radius = dx_ / sqrtf(2.0);
  parallel_for(0, static_cast<int>(nj_), [&](int j) {
    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
      // Estimate from particles
      scalar min_liquid_phi = dx_;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar phi_temp = (pos - p->x_).norm() - std::max(p->radii_, min_radius);
          min_liquid_phi = std::min(min_liquid_phi, phi_temp);
        }
      });

      // "extrapolate" phi into solids if nearby
      scalar solid_phi_val = solid_distance(pos);
      liquid_phi_(i, j) = std::min(min_liquid_phi, solid_phi_val);
    }
  });
}

/*!
  \brief  Add a tracer particle for visualization
*/
void FluidSim::add_particle(const Particle& p) { particles_.push_back(p); }

/*!
  \brief  Particles can still occasionally leave the domain due to truncation
          errors, interpolation error, or large timesteps, so we project them
          back in for good measure.
*/
void FluidSim::particle_boundary_collision(scalar dt) {
  parallel_for(0, static_cast<int>(particles_.size()), [&](int p) {
    Vector2s pp = (particles_[p].x_ - origin_) / dx_;

    // Try commenting this section out to see the degree of accumulated error.
    scalar phi_value = interpolate_value(pp, nodal_solid_phi_);
    if (phi_value < 0) {
      Vector2s normal;
      interpolate_gradient(normal, pp, nodal_solid_phi_);
      normal.normalize();
      particles_[p].x_ -= phi_value * normal;
    }
  });

  m_sorter_->sort(particles_, origin_, dx_);

  auto remover = [&](const Particle& p) {
    return p.x_(0) < origin_(0) - 0.5 * dx_ || p.x_(0) > origin_(0) + (static_cast<scalar>(ni_) + 1.5) * dx_ || p.x_(1) < origin_(1) - 0.5 * dx_ ||
           p.x_(1) > origin_(1) + (static_cast<scalar>(nj_) + 1.5) * dx_;
  };
  particles_.erase(std::remove_if(particles_.begin(), particles_.end(), remover), particles_.end());
}

Vector2s FluidSim::get_velocity_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Vector2s ret = Vector2s::Zero();
  Vector2i ip0 = Vector2i(static_cast<int>(p0(0)), static_cast<int>(p0(1)));
  for (int i = ip0(0) - 1; i <= ip0(0) + 2; ++i) {
    for (int j = ip0(1) - 1; j <= ip0(1) + 2; ++j) {
      if (i < 0 || i > ni_ || j < 0 || j >= nj_) {
        continue;
      }
      Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret(0) += uu(i, j) * w;
    }
  }

  Vector2i ip1 = Vector2i(static_cast<int>(p1(0)), static_cast<int>(p1(1)));
  for (int i = ip1(0) - 1; i <= ip1(0) + 2; ++i) {
    for (int j = ip1(1) - 1; j <= ip1(1) + 2; ++j) {
      if (i < 0 || i >= ni_ || j < 0 || j > nj_) {
        continue;
      }
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret(1) += vv(i, j) * w;
    }
  }

  return ret;
}

Matrix2s FluidSim::get_affine_matrix_quadratic_impl(const Vector2s& position, const Array2s& uu, const Array2s& vv) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s ret = Matrix2s::Zero();
  scalar invD = 4.0 / (dx_ * dx_);
  Vector2i ip0 = Vector2i(static_cast<int>(p0(0)), static_cast<int>(p0(1)));
  for (int i = ip0(0) - 1; i <= ip0(0) + 2; ++i) {
    for (int j = ip0(1) - 1; j <= ip0(1) + 2; ++j) {
      if (i < 0 || i > ni_ || j < 0 || j >= nj_) {
        continue;
      }
      Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret.block<2, 1>(0, 0) += u_(i, j) * w * (position - pos) * invD;
    }
  }

  Vector2i ip1 = Vector2i(static_cast<int>(p1(0)), static_cast<int>(p1(1)));
  for (int i = ip1(0) - 1; i <= ip1(0) + 2; ++i) {
    for (int j = ip1(1) - 1; j <= ip1(1) + 2; ++j) {
      if (i < 0 || i >= ni_ || j < 0 || j > nj_) {
        continue;
      }
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar w = kernel::quadratic_kernel(position - pos, dx_);
      ret.block<2, 1>(0, 1) += v_(i, j) * w * (position - pos) * invD;
    }
  }

  return ret;
}

Vector2s FluidSim::get_velocity_quadratic(const Vector2s& position) { return get_velocity_quadratic_impl(position, u_, v_); }

Matrix2s FluidSim::get_affine_matrix_quadratic(const Vector2s& position) { return get_affine_matrix_quadratic_impl(position, u_, v_); }

Vector2s FluidSim::get_saved_velocity_quadratic(const Vector2s& position) { return get_velocity_quadratic_impl(position, saved_u_, saved_v_); }

Matrix2s FluidSim::get_saved_affine_matrix_quadratic(const Vector2s& position) { return get_affine_matrix_quadratic_impl(position, saved_u_, saved_v_); }

Vector2s FluidSim::get_velocity_and_affine_matrix_with_order(const Vector2s& position, scalar dt, FluidSim::VELOCITY_ORDER v_order,
                                                             FluidSim::INTERPOLATION_ORDER i_order, Matrix2s* affine_matrix) {
  auto get_velocity_func = (i_order == IO_LINEAR) ? &FluidSim::get_velocity : &FluidSim::get_velocity_quadratic;
  auto get_affine_matrix_func = (i_order == IO_LINEAR) ? &FluidSim::get_affine_matrix : &FluidSim::get_affine_matrix_quadratic;

  switch (v_order) {
    case FluidSim::VO_EULER:
      if (affine_matrix) {
        *affine_matrix = (this->*get_affine_matrix_func)(position);
      }
      return (this->*get_velocity_func)(position);
    case FluidSim::VO_RA2: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 2.0 / 3.0 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 2.0 / 3.0 * dt);
        *affine_matrix = a0 * 0.25 + a1 * 0.75;
      }
      return v0 * 0.25 + v1 * 0.75;
    }
    case FluidSim::VO_RK3: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + (v1 * 2.0 - v0) * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + (v1 * 2.0 - v0) * dt);
        *affine_matrix = a0 / 6.0 + a1 * (2.0 / 3.0) + a2 / 6.0;
      }
      return v0 / 6.0 + v1 * (2.0 / 3.0) + v2 / 6.0;
    }
    case FluidSim::VO_RK4: {
      Vector2s v0 = (this->*get_velocity_func)(position);
      Vector2s v1 = (this->*get_velocity_func)(position + v0 * 0.5 * dt);
      Vector2s v2 = (this->*get_velocity_func)(position + v1 * 0.5 * dt);
      Vector2s v3 = (this->*get_velocity_func)(position + v2 * dt);
      if (affine_matrix) {
        Matrix2s a0 = (this->*get_affine_matrix_func)(position);
        Matrix2s a1 = (this->*get_affine_matrix_func)(position + v0 * 0.5 * dt);
        Matrix2s a2 = (this->*get_affine_matrix_func)(position + v1 * 0.5 * dt);
        Matrix2s a3 = (this->*get_affine_matrix_func)(position + v2 * dt);
        *affine_matrix = a0 / 6.0 + a1 / 3.0 + a2 / 3.0 + a3 / 6.0;
      }
      return v0 / 6.0 + v1 / 3.0 + v2 / 3.0 + v3 / 6.0;
    }
    default:
      return Vector2s::Zero();
  }
}

/*!
  \brief  Interpolate velocity from the MAC grid.
*/
Vector2s FluidSim::get_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_ and v_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, u_);
  scalar v_value = interpolate_value(p1, v_);

  return Vector2s(u_value, v_value);
}

Matrix2s FluidSim::get_affine_matrix(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s c_;
  c_.col(0) = affine_interpolate_value(p0, u_) / dx_;
  c_.col(1) = affine_interpolate_value(p1, v_) / dx_;

  return c_;
}

Matrix2s FluidSim::get_saved_affine_matrix(const Vector2s& position) {
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);

  Matrix2s c_;
  c_.col(0) = affine_interpolate_value(p0, saved_u_);
  c_.col(1) = affine_interpolate_value(p1, saved_v_);

  return c_;
}

Vector2s FluidSim::get_saved_velocity(const Vector2s& position) {
  // Interpolate the velocity from the u_ and v_ grids
  Vector2s p = (position - origin_) / dx_;
  Vector2s p0 = p - Vector2s(0, 0.5);
  Vector2s p1 = p - Vector2s(0.5, 0);
  scalar u_value = interpolate_value(p0, saved_u_);
  scalar v_value = interpolate_value(p1, saved_v_);

  return Vector2s(u_value, v_value);
}

Vector2s FluidSim::get_saved_velocity_with_order(const Vector2s& position, FluidSim::INTERPOLATION_ORDER i_order) {
  return (i_order == IO_LINEAR) ? get_saved_velocity(position) : get_saved_velocity_quadratic(position);
}

/*!
  \brief  Given two signed distance values, determine what fraction of a connecting
          segment is "inside"
*/
scalar fraction_inside(scalar phi_left, scalar phi_right) {
  if (phi_left < 0 && phi_right < 0) return 1;
  if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
  if (phi_left >= 0 && phi_right < 0)
    return phi_right / (phi_right - phi_left);
  else
    return 0;
}

/*!
  \brief  Compute finite-volume style face-weights for fluid from nodal signed
          distances
*/
void FluidSim::compute_weights() {
  parallel_for(0, std::max(u_weights_.nj, v_weights_.nj), [this](int j) {
    if (j < u_weights_.nj) {
      for (int i = 0; i < u_weights_.ni; ++i) {
        u_weights_(i, j) = 1 - fraction_inside(nodal_solid_phi_(i, j + 1), nodal_solid_phi_(i, j));
        u_weights_(i, j) = clamp(u_weights_(i, j), 0.0f, 1.0f);
      }
    }
    if (j < v_weights_.nj) {
      for (int i = 0; i < v_weights_.ni; ++i) {
        v_weights_(i, j) = 1 - fraction_inside(nodal_solid_phi_(i + 1, j), nodal_solid_phi_(i, j));
        v_weights_(i, j) = clamp(v_weights_(i, j), 0.0f, 1.0f);
      }
    }
  });
}

/*!
  \brief  An implementation of the variational pressure projection solve for static
          geometry
*/
void FluidSim::solve_pressure(scalar dt) {
  // This linear system could be simplified, but I've left it as is for clarity
  // and consistency with the standard naive discretization
  int system_size = ni_ * nj_;
  if (rhs_.size() != system_size) {
    rhs_.resize(system_size);
    pressure_.resize(system_size);
    matrix_.resize(system_size);
  }
  matrix_.zero();

  // Build the linear system for pressure
  parallel_for(1, nj_ - 1, [&](int j) {
    for (int i = 1; i < ni_ - 1; ++i) {
      int index = i + ni_ * j;
      rhs_[index] = 0;
      pressure_[index] = 0;
      float centre_phi = liquid_phi_(i, j);
      if (centre_phi < 0 && (u_weights_(i, j) > 0.0 || u_weights_(i + 1, j) > 0.0 || v_weights_(i, j) > 0.0 || v_weights_(i, j + 1) > 0.0)) {
        // right neighbour
        float term = u_weights_(i + 1, j) * dt / sqr(dx_);
        float right_phi = liquid_phi_(i + 1, j);
        if (right_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + 1, -term);
        } else {
          float theta = fraction_inside(centre_phi, right_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] -= u_weights_(i + 1, j) * u_(i + 1, j) / dx_;

        // left neighbour
        term = u_weights_(i, j) * dt / sqr(dx_);
        float left_phi = liquid_phi_(i - 1, j);
        if (left_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - 1, -term);
        } else {
          float theta = fraction_inside(centre_phi, left_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += u_weights_(i, j) * u_(i, j) / dx_;

        // top neighbour
        term = v_weights_(i, j + 1) * dt / sqr(dx_);
        float top_phi = liquid_phi_(i, j + 1);
        if (top_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index + ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, top_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] -= v_weights_(i, j + 1) * v_(i, j + 1) / dx_;

        // bottom neighbour
        term = v_weights_(i, j) * dt / sqr(dx_);
        float bot_phi = liquid_phi_(i, j - 1);
        if (bot_phi < 0) {
          matrix_.add_to_element(index, index, term);
          matrix_.add_to_element(index, index - ni_, -term);
        } else {
          float theta = fraction_inside(centre_phi, bot_phi);
          if (theta < 0.01) theta = 0.01;
          matrix_.add_to_element(index, index, term / theta);
        }
        rhs_[index] += v_weights_(i, j) * v_(i, j) / dx_;
      }
    }
  });

  // Solve the system using Robert Bridson's incomplete Cholesky PCG solver
  scalar residual;
  int iterations;
  bool success = solver_.solve(matrix_, rhs_, pressure_, residual, iterations);
  if (!success) {
    std::cout << "WARNING: Pressure solve failed! residual = " << residual << ", iters = " << iterations << std::endl;
  }

  // Apply the velocity update
  parallel_for(0, std::max(u_.nj, v_.nj - 1), [&](int j) {
    if (j < u_.nj) {
      for (int i = 1; i < u_.ni - 1; ++i) {
        int index = i + j * ni_;
        if (u_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i - 1, j) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i - 1, j) >= 0) theta = fraction_inside(liquid_phi_(i - 1, j), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            u_(i, j) -= dt * (pressure_[index] - pressure_[index - 1]) / dx_ / theta;
            u_valid_(i, j) = 1;
          } else {
            u_valid_(i, j) = 0;
          }
        } else {
          u_(i, j) = 0;
          u_valid_(i, j) = 0;
        }
      }    
    }
    if (j >= 1 && j < v_.nj - 1) {
      for (int i = 0; i < v_.ni; ++i) {
        int index = i + j * ni_;
        if (v_weights_(i, j) > 0) {
          if (liquid_phi_(i, j) < 0 || liquid_phi_(i, j - 1) < 0) {
            float theta = 1;
            if (liquid_phi_(i, j) >= 0 || liquid_phi_(i, j - 1) >= 0) theta = fraction_inside(liquid_phi_(i, j - 1), liquid_phi_(i, j));
            if (theta < 0.01) theta = 0.01;
            v_(i, j) -= dt * (pressure_[index] - pressure_[index - ni_]) / dx_ / theta;
            v_valid_(i, j) = 1;
          } else {
            v_valid_(i, j) = 0;
          }
        } else {
          v_(i, j) = 0;
          v_valid_(i, j) = 0;
        }
      }
    }
  });
}

scalar FluidSim::solid_distance(const Vector2s& pos, const Boundary& b) const {
  switch (b.type_) {
    case BT_BOX:
      return b.sign_ * box_distance(pos, b.center_, b.parameter_);
    case BT_CIRCLE:
      return b.sign_ * circle_distance(pos, b.center_, b.parameter_(0));
    case BT_TORUS:
      return b.sign_ * torus_distance(pos, b.center_, b.parameter_(0), b.parameter_(1));
    case BT_TRIANGLE:
      return b.sign_ * triangle_distance(pos, b.center_, b.parameter_(0));
    case BT_HEXAGON:
      return b.sign_ * hexagon_distance(pos, b.center_, b.parameter_(0));
    case BT_CYLINDER:
      return b.sign_ * cylinder_distance(pos, b.center_, b.parameter_(0), b.parameter_(1));
    case BT_UNION:
      return union_distance(solid_distance(pos, *b.op0_), solid_distance(pos, *b.op1_));
    case BT_INTERSECTION:
      return intersection_distance(solid_distance(pos, *b.op0_), solid_distance(pos, *b.op1_));
    default:
      return 1e+20;
  }
}

scalar FluidSim::solid_distance(const Vector2s& pos) const { return solid_distance(pos, *root_boundary_); }

void FluidSim::init_random_particles() {
  int num_particle = ni_ * nj_;
  for (int i = 0; i < ni_; ++i) {
    for (int j = 0; j < nj_; ++j) {
      for (int k = 0; k < 2; ++k) {
        scalar x_ = (static_cast<scalar>(i) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        scalar y = (static_cast<scalar>(j) + 0.5 + ((static_cast<scalar>(rand()) / static_cast<scalar>(RAND_MAX)) * 2.0 - 1.0)) * dx_;
        Vector2s pt = Vector2s(x_, y) + origin_;

        scalar phi = solid_distance(pt);
        if (phi > dx_ * 20.0) particles_.emplace_back(pt, Vector2s::Zero(), dx_ / sqrt(2.0), rho_);
      }
    }
  }
}

void FluidSim::map_p2g() {
  if (interpolation_order == IO_LINEAR)
    map_p2g_linear();
  else
    map_p2g_quadratic();
}

void FluidSim::map_p2g_linear() {
  // u-component of velocity
  parallel_for(0, nj_ + 1, [this](int j) {
    if (j < nj_) {
      for (int i = 0; i < ni_ + 1; ++i) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        scalar sumw = 0.0;
        scalar sumu = 0.0;
        m_sorter_->getNeigboringParticles_cell(i, j, -1, 0, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
            sumu += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
            sumw += w;
          }
        });

        u_(i, j) = sumw ? sumu / sumw : 0.0;
      }
    }

    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar sumw = 0.0;
      scalar sumu = 0.0;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -1, 0, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar w = p->mass_ * kernel::linear_kernel(p->x_ - pos, dx_);
          sumu += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
          sumw += w;
        }
      });

      v_(i, j) = sumw ? sumu / sumw : 0.0;
    }
  });
}

void FluidSim::map_p2g_quadratic() {
  // u-component of velocity
  parallel_for(0, nj_ + 1, [this](int j) {
    if (j < nj_) {
      for (int i = 0; i < ni_ + 1; ++i) {
        Vector2s pos = Vector2s(i * dx_, (j + 0.5) * dx_) + origin_;
        scalar sumw = 0.0;
        scalar sumu = 0.0;
        m_sorter_->getNeigboringParticles_cell(i, j, -2, 1, -1, 1, [&](const NeighborParticlesType& neighbors) {
          for (const Particle* p : neighbors) {
            scalar w = p->mass_ * kernel::quadratic_kernel(p->x_ - pos, dx_);
            sumu += w * (p->v_(0) + p->c_.col(0).dot(pos - p->x_));
            sumw += w;
          }
        });

        u_(i, j) = sumw > 0.0 ? sumu / sumw : 0.0;
      }
    }

    for (int i = 0; i < ni_; ++i) {
      Vector2s pos = Vector2s((i + 0.5) * dx_, j * dx_) + origin_;
      scalar sumw = 0.0;
      scalar sumu = 0.0;
      m_sorter_->getNeigboringParticles_cell(i, j, -1, 1, -2, 1, [&](const NeighborParticlesType& neighbors) {
        for (const Particle* p : neighbors) {
          scalar w = p->mass_ * kernel::quadratic_kernel(p->x_ - pos, dx_);
          sumu += w * (p->v_(1) + p->c_.col(1).dot(pos - p->x_));
          sumw += w;
        }
      });

      v_(i, j) = sumw > 0.0 ? sumu / sumw : 0.0;
    }
  });
}

/*!
  \brief  A general affine FLIP scheme that unifies all the other FLIP schemes
          used in this code
*/
void FluidSim::map_g2p_flip_general(float dt, const scalar lagrangian_ratio, const scalar lagrangian_symplecticity, const scalar affine_stretching_ratio,
                                    const scalar affine_rotational_ratio) {
  bool use_affine = affine_stretching_ratio > 0. || affine_rotational_ratio > 0.;
  parallel_for(0, static_cast<int>(particles_.size()), [&](int i) {
    auto& p = particles_[i];

    Matrix2s C = Matrix2s::Zero();
    Vector2s next_grid_velocity = get_velocity_and_affine_matrix_with_order(p.x_, dt, velocity_order, interpolation_order, use_affine ? (&C) : nullptr);
    Vector2s lagrangian_velocity = p.v_;
    Vector2s original_grid_velocity;

    if (lagrangian_ratio > 0.0) {
      original_grid_velocity = get_saved_velocity_with_order(p.x_, interpolation_order);
      p.v_ = next_grid_velocity + (lagrangian_velocity - original_grid_velocity) * lagrangian_ratio;
    } else {
      p.v_ = next_grid_velocity;
    }

    p.c_ = C * (affine_stretching_ratio + affine_rotational_ratio) * 0.5 + C.transpose() * (affine_stretching_ratio - affine_rotational_ratio) * 0.5;

    if (lagrangian_symplecticity > 0.0 && lagrangian_ratio > 0.0) {
      scalar logdJ = C.trace() * dt;
      p.logJ_ += logdJ;
      scalar pos_adj = lagrangian_symplecticity;
      Vector2s pp = (p.x_ - origin_) / dx_;
      scalar phi_value = interpolate_value(pp, nodal_solid_phi_);
      if (phi_value < 0) {
        Vector2s normal;
        interpolate_gradient(normal, pp, nodal_solid_phi_);
        normal.normalize();
        float dotnv = -normal.dot(lagrangian_velocity);
        if (dotnv > 0.0) {
          pos_adj = 0.0;
        }
      }
      if (p.logJ_ < -0.001) {
        pos_adj = 0.0;
      }
      p.x_ += (next_grid_velocity + (lagrangian_velocity - original_grid_velocity) * lagrangian_ratio * pos_adj) * dt;
    } else {
      p.x_ += next_grid_velocity * dt;
    }
  });
}

void FluidSim::render_boundaries(const Boundary& b) {
  switch (b.type_) {
    case BT_CIRCLE:
      draw_circle2d(b.center_, b.parameter_(0), 64);
      break;
    case BT_BOX:
      draw_box2d(b.center_, b.parameter_(0), b.parameter_(1));
      break;
    case BT_TORUS:
      draw_circle2d(b.center_, b.parameter_(0), 64);
      draw_circle2d(b.center_, b.parameter_(1), 64);
      break;
    case BT_HEXAGON:
      draw_circle2d(b.center_, b.parameter_(0), 6);
      break;
    case BT_TRIANGLE:
      draw_circle2d(b.center_, b.parameter_(0), 3);
      break;
    case BT_UNION:
    case BT_INTERSECTION:
      render_boundaries(*b.op0_);
      render_boundaries(*b.op1_);
      break;
    default:
      break;
  }
}

void FluidSim::render() {
  glPushMatrix();
  glScaled(1.0 / (dx_ * ni_), 1.0 / (dx_ * ni_), 1.0 / (dx_ * ni_));

  if (draw_grid_) {
    glColor3f(0.75, 0.75, 0.75);
    glLineWidth(1);
    draw_grid2d(origin_, dx_, ni_, nj_);
  }

  if (draw_boundaries_) {
    render_boundaries(*root_boundary_);
  }

  if (draw_velocities_) {
    glColor3f(1, 0, 0);
    scalar crit = dx_ * dx_ * 100.0f;
    glBegin(GL_LINES);
    for (int j = 0; j < nj_; ++j) {
      for (int i = 0; i < ni_; ++i) {
        Vector2s pos = Vector2s((i + 0.5) * dx_, (j + 0.5) * dx_) + origin_;
        Vector2s vel = get_velocity(pos);
        if (vel.squaredNorm() > crit) {
          Vector2s pos_v = pos + 0.01 * vel;
          glVertex2fv(pos.data());
          glVertex2fv(pos_v.data());
        }
      }
    }
    glEnd();
  }

  if (draw_particles_) {
    glColor3f(0, 0, 1);
    glPointSize(5);
    glBegin(GL_POINTS);
    for (unsigned int i = 0; i < particles_.size(); ++i) {
      glVertex2fv(particles_[i].x_.data());
    }
    glEnd();
  }

  glPopMatrix();
}

FluidSim::Boundary::Boundary(const Vector2s& center, const Vector2s& parameter, BOUNDARY_TYPE type, bool inside)
    : center_(center), parameter_(parameter), type_(type), sign_(inside ? -1.0 : 1.0) {}

FluidSim::Boundary::Boundary(Boundary* op0, Boundary* op1, BOUNDARY_TYPE type) : op0_(op0), op1_(op1), type_(type), sign_(op0 ? op0->sign_ : false) {}

Particle::Particle(const Vector2s& x, const Vector2s& v, const scalar& radii, const scalar& density)
    : x_(x), v_(v), radii_(radii), mass_(4.0 / 3.0 * M_PI * radii * radii * radii * density), logJ_(0) {
  c_.setZero();
  buf0_.setZero();
}

Particle::Particle() : x_(Vector2s::Zero()), v_(Vector2s::Zero()), radii_(0.0), mass_(0.0), logJ_(0.0) {
  c_.setZero();
  buf0_.setZero();
}

Particle::Particle(const Particle& p) : x_(p.x_), v_(p.v_), radii_(p.radii_), mass_(p.mass_), logJ_(0.0) {
  c_.setZero();
  buf0_.setZero();
}

/*!
  \brief Apply several iterations of a very simple "Jacobi"-style propagation of
         valid velocity data in all directions
*/
void extrapolate(Array2s& grid, Array2s& old_grid, const Array2s& grid_weight, const Array2s& grid_liquid_weight, Array2c& valid, Array2c old_valid_,
                 const Vector2i& offset, int num_layers) {
  if (!num_layers) return;

  // Initialize the list of valid_ cells
  for (int j = 0; j < valid.nj; ++j) valid(0, j) = valid(valid.ni - 1, j) = 0;
  for (int i = 0; i < valid.ni; ++i) valid(i, 0) = valid(i, valid.nj - 1) = 0;

  Array2s* pgrids[] = {&grid, &old_grid};
  Array2c* pvalids[] = {&valid, &old_valid_};

  for (int j = 1; j < grid.nj - 1; ++j) {
    for (int i = 1; i < grid.ni - 1; ++i) {
      valid(i, j) = grid_weight(i, j) > 0 && (grid_liquid_weight(i, j) < 0 || grid_liquid_weight(i + offset(0), j + offset(1)) < 0);
    }
  }

  old_valid_ = valid;
  old_grid = grid;

  for (int layers = 0; layers < num_layers; ++layers) {
    Array2s* pgrid_source = pgrids[layers & 1];
    Array2s* pgrid_target = pgrids[!(layers & 1)];

    Array2c* pvalid_source = pvalids[layers & 1];
    Array2c* pvalid_target = pvalids[!(layers & 1)];

    parallel_for(1, grid.nj - 1, [&](int j) {
      for (int i = 1; i < grid.ni - 1; ++i) {
        scalar sum = 0;
        int count = 0;

        if (!(*pvalid_source)(i, j)) {
          if ((*pvalid_source)(i + 1, j)) {
            sum += (*pgrid_source)(i + 1, j);
            ++count;
          }
          if ((*pvalid_source)(i - 1, j)) {
            sum += (*pgrid_source)(i - 1, j);
            ++count;
          }
          if ((*pvalid_source)(i, j + 1)) {
            sum += (*pgrid_source)(i, j + 1);
            ++count;
          }
          if ((*pvalid_source)(i, j - 1)) {
            sum += (*pgrid_source)(i, j - 1);
            ++count;
          }

          // If any of neighbour cells were valid_,
          // assign the cell their average value and tag it as valid_
          if (count > 0) {
            (*pgrid_target)(i, j) = sum / static_cast<scalar>(count);
            (*pvalid_target)(i, j) = 1;
          }
        }
      }
    });

    if (layers != num_layers - 1) {
      for (int j = 1; j < grid.nj - 1; ++j) {
        for (int i = 1; i < grid.ni - 1; ++i) {
          if (!(*pvalid_source)(i, j)) {
            (*pgrid_source)(i, j) = (*pgrid_target)(i, j);
            (*pvalid_source)(i, j) = (*pvalid_target)(i, j);
          }
        }
      }
    }
  }
}
