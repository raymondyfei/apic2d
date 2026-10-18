Affine Particle in Cell in 2D
================
![apic2d_pad](https://github.com/nepluno/apic2d/assets/5384915/e0bf3f6b-85ef-4809-9911-c786800ea79f)

APIC2D is an educational project to illustrate the affine-particle-in-cell algorithm in 2D for water simulation. The code can be compiled and run on Windows, Mac OS X (either Intel or Apple Silicon), Linux, and [webpage](https://nepluno.github.io/apic2d/) (through [WASM, or WebAssembly](https://webassembly.org/)). It is well-optimized. For the default setting, the simulation runs in real-time on an iPad mini 6.

The papers implemented here include:

Jiang, Chenfanfu, et al. "The affine particle-in-cell method." ACM Transactions on Graphics (TOG) 34.4 (2015): 51.

Batty, Christopher, Florence Bertails, and Robert Bridson. "A fast variational framework for accurate solid-fluid coupling." ACM Transactions on Graphics (TOG). Vol. 26. No. 3. ACM, 2007.

Ando, Ryoichi, Nils Thurey, and Reiji Tsuruno. "Preserving fluid sheets with adaptively sampled anisotropic particles." IEEE transactions on visualization and computer graphics 18.8 (2012): 1202-1214.

Brackbill, Jeremiah U., and Hans M. Ruppel. "FLIP: A method for adaptively zoned, particle-in-cell calculations of fluid flows in two dimensions." Journal of Computational physics 65.2 (1986): 314-343.

Fei, Yun (Raymond), et al. "Revisiting Integration in the Material Point Method: A Scheme for Easier Separation and Less Dissipation." ACM Transactions on Graphics (TOG) 40.4 (2021): 109.

It contains multiple integrators that you may switch and compare through changing the `integration_scheme` variable in the code. Its value can be one of the following:
```
IT_PIC: original particle-in-cell (PIC)
IT_FLIP: original fluid-implicit-particle (FLIP)
IT_RPIC: rotational particle-in-cell (RPIC)
IT_APIC: affine particle-in-cell (APIC)
IT_AFLIP: affine fluid-implicit-particle (AFLIP)
IT_ASFLIP: affine separable fluid-implicit-particle (ASFLIP)
```

It also supports using different orders for velocity evaluation, where one may change the `velocity_order` variable in the code. Its value can be one of the following:
```
VO_EULER: first order evaluation
VO_RA2: Ralston's second order evaluation
VO_RK3: Runge Kutta's 3rd-order method
VO_RK4: Runge Kutta's 4rd-order method
```

Dependencies
--------------------
APIC2D depends on the Eigen libraries (included), TBB (optional but provides the best performance), and GLUT/FreeGLUT for simple visualization. 

For Windows or WASM, CMake should automatically handle the dependencies. For Mac OS X, TBB needs to be installed, which can be done with Homebrew. For Ubuntu Linux, you may use the APT package handling utility to install these dependencies. 

Compilation
-----------------
To compile APIC2D, you'll need CMake or CMake-GUI (https://cmake.org). Additionally, for WASM, please follow the documentation of [Emscripten SDK](https://emscripten.org/docs/getting_started/index.html).

CMake:
1. make a directory, say, `build`, with `mkdir build`
2. enter the `build` directory. For Desktop, type `cmake ..` (or `cmake -G Xcode ..` to generate Xcode project files on Mac; for Windows please use `cmake -G <generator> ..` to use the specific generators, or simply type `cmake -G` to list all the available generators). For WASM type `emcmake cmake ..`.
3. Optionally you can adjust the options with `ccmake .`
4. type `make` to compile the code. You may use `make -j` to speed up the compilation process.

CMake-GUI:
1. open CMake-GUI, enter the correct directory for source code, and build. Then click `Configure`, choose the generator (for Windows, select the installed version of the Microsoft Visual Studio).
2. CMake should automatically find all the dependencies. If not, check the `Advanced` box and locate those missing libraries manually. On Windows, please make sure you have picked the libraries corresponding to the architecture you have selected (say, 32-bit libraries for x86 and 64-bit libraries for x64).
3. click generate after fixing all missing variables.
4. open the solution (for Visual Studio or Xcode) and compile the code.

User Interface
-----------------
Run `./build/apic2d/APIC2D --timing` to print average phase timings every 100
simulation steps. The `[pcg iterate]` entry includes total iterations,
iterations per solve and per step, average milliseconds per iteration, and a
breakdown of matrix multiplication, preconditioner application, and vector
operations/reductions in milliseconds per step. These nested timings describe
the iteration loop only (not the initial preconditioner application) and are
not added again to `[total]`. Detailed profiling adds timer overhead and is
disabled without `--timing`.

Press `<LMB>` (or finger touch on a mobile device) and drag directly on the simulated region to paint velocities.

Press `<shift>` with mouse buttons to move the camera.
1. Pan: `<shift>+<LMB>`
2. Zoom-in (drag to mark the region): `<shift>+<MMB>`
3. Zoom-out: `<shift>+<RMB>`

Optimizations
-----------------
The following are approximate local measurements from incremental comparisons,
not guaranteed speedups. Each percentage describes the stated metric relative
to that experiment's baseline; the percentages are **not additive**. Results
depend on fluid occupancy, resolution, thread count, and hardware.
The list is ordered roughly by practical impact on overall runtime, not by the
largest percentage reduction in an individual kernel; the different benchmark
baselines prevent a strict like-for-like ranking.

| Optimization | Measured reduction | Benchmark context |
|---|---|---|
| Active-fluid-cell compaction | **~40% less pressure-phase time; ~20% less whole-step time** | 100x100 grid, 300 steps, one CPU thread, ~12-14% active cells. Includes mapping and scattering overhead. |
| Specialized MIC application | **~40% less MIC-application time; ~14-17% less pressure-phase time** | 100x100 grid, 300 steps, one CPU thread. Includes layout conversion; both versions took 4,172 total iterations. |
| Frontier-based velocity extrapolation | **~32% less extrapolation time** | Latest combined 100x100, four-thread run: 0.0830 to 0.0562 ms/step. |
| Unrolled fixed-width matrix multiplication | **~49% less matrix-multiplication time** | Latest combined 100x100, four-thread run: 0.0322 to 0.0163 ms/step. |
| Fused vector updates and SIMD residual norm | **~50% less vector-operation time; ~5% less solver time** | Fixed 100x100 pressure system, 24 iterations: vector work ~0.21 to ~0.10 ms/solve; solver ~1.68 to ~1.59 ms/solve. |
| Serial pressure scattering on small grids | **~69% less pressure-scatter time** | Latest combined 100x100, four-thread run: 0.0124 to 0.0038 ms/step. |
| Per-row assembly reuse and direct stencil-slot writes | **~30% less assembly time** | Moving-fluid assembly ~0.054 to ~0.038 ms/step. Whole-step improvement was too small to distinguish confidently from noise. |
| Pressure warm-starting | **~4% fewer PCG iterations** | 100x100 grid, 300 steps: 14.40 to 13.84 iterations/solve. This is not a wall-time reduction; initial-residual computation adds work. |

The fixed-width matrix, serial-scatter, and frontier-extrapolation rows are phase
measurements from the same combined comparison, not isolated whole-step gains
for each feature. Together, these changes reduced
whole-step time by **~11%** on the 100x100, four-thread, 300-step sample
(0.478 to 0.427 ms/step), and **~8%** on a 256x256, four-thread, 100-step sample
(8.97 to 8.27 ms/step). The 100x100 single-thread sample improved by **~2%**
(0.886 to 0.866 ms/step). Floating-point ordering differences can slightly change
iteration counts and subsequent fluid motion.

Whole-matrix CSR/MIC structure caching was tried and removed: the moving-fluid
benchmark had **0% whole-matrix reuse hits**. Per-row assembly reuse remains.
Serial assembly and velocity updates were also tried, but were slower and
remain parallel. Profiling and residual-norm correctness fixes are not counted
as speedups.

### Implementation details

Pressure solves warm-start from the previous pressure by default, clearing
inactive cells and falling back to a zero guess if the initial residual is
larger than the right-hand-side norm. The convergence tolerance remains
`1e-4 * max(abs(rhs))`. To compare against zero-initialized solves, run
`./build/apic2d/APIC2D --timing --cold-start`. Compare iteration counts and
overall timings, including the new `[pcg initial residual]` phase, which
accounts for computing and checking the initial guess.

PCG fuses pressure and residual updates with the residual-norm reduction in
one traversal. Standalone residual norms use Eigen's value-only SIMD reduction.
The preconditioner and RHS-relative convergence tolerance are unchanged.

The pressure system contains only active interior fluid cells, in row-major
order. Warm-start values are gathered from the full pressure grid before each
solve and scattered back afterward, with inactive cells set to zero. This
preserves the warm start when the active set changes. With `--timing`, the
report includes average active unknowns versus total grid cells, `[pressure map]`
and `[pressure scatter]` timings, and an inclusive `[pressure solve total]`
covering mapping, assembly, the solver, scattering, and velocity updates. The
inclusive total is a separate summary and is not added again to `[total]`.

Pressure solves apply MIC using fixed-width triangular rows for the compact
five-point stencil, avoiding generic sparse traversal and the forward solve's
RHS copy. The original MIC factorization is retained. Conversion to the
fixed-width layout is included in `[pcg factorize]`. The general-purpose
`PCGSolver` keeps its generic application by default; callers can opt in with
`set_use_grid_preconditioner(true)`. Factors with more than two entries in any
triangular row or column use the generic implementation.

Pressure matrix multiplication also uses unrolled fixed-width rows with up to
five entries. Matrix setup is reported as `[pcg matrix build]`; eligible solves
do not build CSR as well. General-purpose `PCGSolver` callers can opt in with
`set_use_grid_matrix(true)`; matrices with wider rows fall back to CSR.

Pressure assembly retains matching column-index rows, including when other
rows change, and writes coefficients directly into their stencil slots.
Topology changes and fluid-face weights crossing zero rebuild affected rows.
The timing report's `assembly row reuse` line shows the percentage of rows
retained. The selected matrix layout and MIC structures are rebuilt for each
solve that needs them; whole-matrix structure caching is not used.

Pressure scattering runs serially below 65,536 grid cells to avoid scheduling
overhead, and in parallel for larger grids. Assembly and velocity updates stay
parallel. Velocity extrapolation keeps the same two layers of neighbor averaging
but processes a frontier of newly reachable faces instead of repeatedly scanning
and copying entire grids. Each layer is computed before its values become valid
for the next layer. Frontier storage is reused between velocity components and
steps, with the validity mask reset to each component's dimensions.

### Regression checks

The solver regression checks can be run without the visualization dependencies:
```
c++ -std=c++14 -O3 -DSINGLE_THREADED -Iinclude/eigen -Iapic2d tests/pcg_solver_test.cpp -o build/pcg_solver_test
./build/pcg_solver_test
```

On macOS, the compact-projection and extrapolation checks can also be built and run with:
```
c++ -std=c++14 -O3 -DSINGLE_THREADED -Wno-deprecated-declarations -Iinclude/eigen -Iapic2d \
    tests/pressure_projection_test.cpp apic2d/{fluidsim,kernel,sorter,openglutils,gluvi}.cpp \
    -framework OpenGL -framework GLUT -o build/pressure_projection_test
./build/pressure_projection_test
```
