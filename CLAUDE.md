# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System and Commands

### Environment Setup
```bash
# Required: Set Castro path (parent hydrodynamics code)
export CASTRO_HOME=../../../Castro

# Required for implicit radiation solvers (MGFLD, SGFLD)
export HYPRE_DIR=/path/to/hypre
# Or for OpenMP builds:
export HYPRE_OMP_DIR=/path/to/hypre-with-openmp
```

### Building a Test Problem
```bash
cd Exec/RadThermalWave  # or any test problem directory

# Edit GNUmakefile to set:
# - DIM = 1, 2, or 3 (problem dimension)
# - DEBUG = TRUE/FALSE
# - COMP/FCOMP (compiler choices)
# - EOS_dir, Network_dir, Opacity_dir (physics modules)

# Build
gmake -j4

# Executable will be named like: Castro2d.Linux.g++.gfortran.ex
```

### Running Simulations
```bash
# Run with input file and probin file
./Castro2d.Linux.g++.gfortran.ex inputs.2d

# For MPI parallel runs
mpiexec -n 4 ./Castro2d.Linux.g++.gfortran.MPI.ex inputs.2d

# Common runtime parameters to adjust in input files:
# - amr.n_cell: base grid resolution
# - amr.max_level: AMR refinement levels
# - max_step: maximum timesteps
# - amr.plot_int: output frequency
```

### Building Analysis Tools
```bash
# Some test problems have analysis tools
cd Exec/RadSphere/Tools
make programs=fradsphere  # builds fradsphere.exe for post-processing
```

## Architecture Overview

### Core Radiation-Hydrodynamics Coupling
The code implements radiation hydrodynamics with tight coupling between Castro (hydrodynamics) and radiation transport:

1. **Main Driver**: `Radiation` class (Source/Radiation.H/cpp) orchestrates radiation solvers
2. **Solver Types** (set via runtime parameters):
   - `SingleGroupSolver`: Gray (frequency-integrated) radiation
   - `SGFLDSolver`: Single-group flux-limited diffusion
   - `MGFLDSolver`: Multi-group flux-limited diffusion (most sophisticated)
3. **Frame Treatment**: Supports both lab frame and comoving frame radiation

### Modular Source Structure
- **Source/**: Core radiation algorithms (dimension-independent)
  - `RadSolve.cpp`: Main radiation implicit solver interface
  - `MGFLD.cpp/MGFLDRadSolver.cpp`: Multi-group FLD implementation
  - `rad_params.f90`: Physical constants and group structure
- **Source/Src_1d,2d,3d/**: Dimension-specific implementations
  - Separate implementations for efficiency (no runtime dimension checks)
  - Each contains specialized hydro coupling, boundary conditions, flux calculations
- **Exec/*/**: Individual test problems with:
  - `Prob_Nd.f90`: Problem initialization
  - `probdata.f90`: Problem-specific parameters
  - Local `Make.package`: Additional problem-specific sources

### Key Solver Components
1. **Implicit System**: Uses Hypre library for solving radiation diffusion equations
   - `HypreMultiABec`: Multi-group linear system solver
   - Critical for stiff radiation-matter coupling
2. **Flux Limiters** (`fluxlimiter.f90`): Ensures causality in radiation transport
3. **Opacity Models** (`Opacity/*/`): Modular opacity implementations
4. **EOS Integration** (`EOS/*/`): Couples to Castro's equation of state framework

## Test Problems and Verification

### Core Verification Tests
- **RadThermalWave**: Marshak wave problem - tests radiation diffusion
- **RadSphere**: Spherical radiation source - has analytic solution
- **Rad2Tshock**: Radiation hydrodynamic shocks - tests coupling
- **RadSuOlson**: Su-Olson benchmark - standard radiation transport test
- **RadBreakout**: Shock breakout from stellar surfaces

### Running Convergence Studies
```bash
# Typical convergence test workflow
cd Exec/RadSphere
# Run with increasing resolution
for ncell in 32 64 128 256; do
    sed -i "s/amr.n_cell.*/amr.n_cell = $ncell/" inputs
    ./Castro1d.ex inputs
    mv plt* results_${ncell}/
done
# Analyze with python scripts in test_problem/python/
```

## Development Patterns

### Adding New Opacity Models
1. Create directory in `Opacity/your_model/`
2. Implement `opacity_table_module.f90` with required interfaces
3. Set `Opacity_dir := your_model` in GNUmakefile

### Modifying Radiation Solvers
- Changes to algorithms typically go in `Source/RadSolve.cpp` or `Source/MGFLD*.cpp`
- Dimension-specific optimizations belong in `Source/Src_Nd/`
- Boundary conditions: modify `RadBndry.cpp` and dimension-specific `RadBndry_Nd.f90`

### Input File Hierarchy
1. **GNUmakefile**: Build-time physics choices (EOS, network, opacity)
2. **inputs**: Runtime AMR, I/O, and solver parameters
3. **probin**: Problem-specific Fortran namelist parameters

### Debugging Radiation Issues
```bash
# Enable debug build
sed -i 's/DEBUG = FALSE/DEBUG = TRUE/' GNUmakefile
gmake clean; gmake

# Useful runtime parameters for debugging:
# castro.v = 2  # verbose output
# amr.v = 2     # AMR verbosity
# radiation.v = 2  # radiation solver details
# radiation.plot_lambda = 1  # output flux limiter
# radiation.plot_kappa_r = 1  # output opacities
```

## Important Implementation Notes

### Radiation Energy Groups
- Group structure defined at runtime via input parameters
- Neutrino transport uses species-specific groups (electron/anti-electron/heavy)
- Group boundaries typically logarithmic in frequency space

### Boundary Conditions
- Radiation BCs separate from hydro BCs
- Marshak (equilibrium) and vacuum BCs most common
- Set via `radiation.lo_bc` and `radiation.hi_bc` in input files

### Performance Considerations
- MGFLD solver is expensive - use geometric multigrid (`hmabec.verbose = 2` to monitor)
- Comoving frame (`comoving = 1`) more stable but expensive
- Implicit updates dominate cost - tune `radiation.maxiter` carefully

### Common Pitfall Avoidance
- Always ensure `HYPRE_DIR` is set when using implicit solvers
- Check that opacity and EOS models are consistent
- Verify units consistency between radiation and hydro (cgs vs code units)
- For multigroup: ensure group structure captures relevant physics

## Opacity Calculations

The code uses a multigroup radiation transport approach where opacity calculations are central to radiation-matter coupling.

### Multigroup Structure
Each radiation group has:
- `nugroup(g)`: Group center frequency
- `xnu(g), xnu(g+1)`: Group boundary frequencies  
- `dnugroup(g)`: Group frequency width

### Planck Mean Opacity (κ_P)
Used for emission and absorption:
```
κ_P = ∫ κ(ν)B(ν,T)dν / ∫ B(ν,T)dν
```

**Implementation in `ca_compute_emissivity`**:
```fortran
jg(i,j,g) = Bg * kap(i,j,g)  ! Emission coefficient
```
Where:
- `jg`: Group emission coefficient
- `Bg`: Integrated Planck function for group
- `kap`: Group Planck mean opacity

**Planck Function Integration Methods**:
1. **Exact integration** (`integrate_Planck = 1`):
   ```fortran
   call BdBdTIndefInteg(T, xnu(g), B0, dBdT0)
   call BdBdTIndefInteg(T, xnu(g+1), B1, dBdT1)
   Bg = B1 - B0  ! Integrated B for group g
   ```

2. **Group center evaluation** (default):
   ```fortran
   nu = nugroup(g)
   Bg = 8π*h/c³ * nu³/(exp(hν/kT) - 1) * dnu
   ```

3. **Gray radiation** (`ngroups = 1`):
   ```fortran
   Bg = arad * T⁴  ! Stefan-Boltzmann
   ```

### Rosseland Mean Opacity (κ_R)
Used for radiation diffusion:
```
1/κ_R = ∫ (1/κ(ν)) ∂B/∂T dν / ∫ ∂B/∂T dν
```
- Stored as `kpr(i,j,g)` for each group
- Used in flux limiter calculations
- Harmonic mean suitable for diffusion

### Opacity Models

**Simple Power Law** (`powerlaw-example/`):
```fortran
kp = const_kappa_p * ρ^m * T^(-n) * ν^p
kr = const_kappa_r * ρ^m * T^(-n) * ν^p
```

**Thomson Scattering** (`breakout/`):
```fortran
kp = ρYe * 0.4 * fac  ! Planck mean
kr = ρYe * 0.4        ! Rosseland mean (cm²/g)
```

### Key Implementation Files
- `Source/Src_Nd/MGFLD_Nd.f90`: Multigroup opacity calculations
- `Source/blackbody.f90`: Planck function integration
- `Opacity/*/opacity_table_module.f90`: Opacity model implementations
- `Source/Src_Nd/ca_opacs`: Calls to opacity modules

### Temperature Derivatives
For implicit solver stability:
```fortran
call get_opacities(kp1, kr1, ρ, T-dT, Ye, ν, ...)
call get_opacities(kp2, kr2, ρ, T+dT, Ye, ν, ...)
dkdT(i,j,g) = (kp2 - kp1)/(2*dT)
```

### Practical Notes
- Each group calls `get_opacities` with its representative frequency
- Emission uses Kirchhoff's law: `jg = κ_P * B`
- Opacity modules can be swapped via `Opacity_dir` in GNUmakefile
- For debugging: set `radiation.plot_kappa_r = 1` to output opacities