# BeamletOpticsPrescriptionIO.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Julia](https://img.shields.io/badge/Julia-v1.10+-purple.svg)](https://julialang.org)
[![Tested on](https://img.shields.io/badge/Tested%20Models-3%2C700%2B-brightgreen.svg)](#benchmarking--validation)

> **Optical Prescription File Importer (Zemax `.zmx`), Component Catalog Interface, and Standalone Code Generator for [BeamletOptics.jl](https://github.com/7BitLogic/BeamletOptics.jl)**

---

## Overview

**BeamletOpticsPrescriptionIO.jl** bridges digital optical prescription designs (starting with Zemax OpticStudio `.zmx` files) with modern 3D physical ray and beamlet tracing in Julia. It parses sequential surface prescriptions and transforms them into native, volumetric 3D Signed Distance Function (SDF) optical objects for `BeamletOptics.jl`.

### Key Capabilities

- **Robust Zemax Parsing**: Handles UTF-16LE (with or without BOM), UTF-8, and ISO-8859-1 encodings. Automatically extracts curvatures, thicknesses, clear apertures, conic constants, even aspheric polynomial coefficients, and multi-configuration parameters (`THIC`).
- **Sequential-to-Element Assembly**: Groups adjacent sequential optical surfaces into concrete 3D optical assemblies:
  - Singlet lenses (`Lens` with spherical or `EvenAsphericalSurface` profiles)
  - Cemented doublets (`DoubletLens` / `SphericalDoubletLens`)
  - Cemented triplets (`TripletLens` / `SphericalTripletLens`)
  - Reflective mirrors (`RoundPlanoMirror`, `SphericalMirror`)
  - Aperture stops (`ZmxStop`)
  - Image detectors (`Detector`)
- **Direct Catalog Integration (`load_lens_from_zmx_cat`)**:
  - Load individual lenses directly from manufacturer catalogs (e.g., Thorlabs) as native `BeamletOptics` components.
  - Case-insensitive search, automatic 3D axial positioning, and metadata extraction.
- **Glass & Dispersion Engine**:
  - Connects to `RefractiveIndex.jl` for over 1,000 commercial optical glasses (Schott, Ohara, Hoya, CDGM, Sumita).
  - Built-in physical Cauchy dispersion calculation for model glasses with specified Abbe numbers ($V_d$).
  - Extensive historical glass database covering classic patent formulations.
- **Analytical Refocusing & Quality Metrics**:
  - Closed-form least-squares focal plane optimization (`find_best_focus`, `refocus!`).
  - Automated Airy disk radius calculation ($r_{\text{Airy}} = 0.61 \frac{\lambda}{\text{NA}}$) and diffraction-limit verification.
- **Code Generation & CLI**:
  - Generates standalone, human-readable Julia scripts (`.jl`) that can be executed independently without any dependency on the prescription importer.
- **Backwards Compatibility**:
  - Exports `BeamletOpticsZMX` as an alias for backwards compatibility with earlier scripts and workflows.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/7BitLogic/BeamletOpticsPrescriptionIO.jl")
```

Or for local development:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

---

## Minimal Working Examples (MWE)

### 1. Load an Individual Lens from a Catalog

Load a commercial achromatic doublet or singlet directly by part number and place it along the optical axis ($+Y$ in BeamletOptics):

```julia
using BeamletOptics
using BeamletOpticsPrescriptionIO

# Load a Thorlabs 1-inch, 50 mm achromatic doublet and position it at y = 100 mm
achromat = load_lens_from_zmx_cat(:thorlabs, "AC254-050-A", position=0.10)

# Build a simple optical system with a detector in the focal plane (EFL ≈ 50 mm)
detector = Detector(0.025)
translate3d!(detector, [0.0, 0.10 + 0.0432, 0.0]) # Back focal length ~43.2 mm

sys = System([achromat, detector])

# Trace a collimated beam (15 mm diameter, 587.6 nm d-line)
src = CollimatedSource([0.0, 0.05, 0.0], [0.0, 1.0, 0.0], 0.015, 587.56e-9; num_rays=80, num_rings=4)
solve_system!(sys, src)

hits = spot_diagram(detector)
println("Detected ray hits: ", length(hits))
```

### 2. Import a Complete Zemax System

Import an entire multi-element optical prescription (e.g. Double Gauss camera lens) and trace rays:

```julia
using BeamletOptics
using BeamletOpticsPrescriptionIO

# Import full system prescription
res = import_zmx("path/to/lens_design.zmx")

# Inspect parsed elements
println("Optical elements found: ", length(res.elements))

# Automatically configure ray source matching the entrance pupil diameter and conjugate
source = suggest_source(res; num_rays=200, num_rings=6)

# Trace through the system
solve_system!(res.system, source)

# Evaluate spot diagram
hits = spot_diagram(res.detector)
println("Hits on image sensor: ", length(hits))
```

### 3. Automated Focus Optimization & Refocusing

Analytically determine the optimal focal plane position to minimize the geometric RMS spot radius:

```julia
using BeamletOptics
using BeamletOpticsPrescriptionIO

res = import_zmx("path/to/lens_design.zmx")
source = suggest_source(res)
solve_system!(res.system, source)

# Find optimal focus displacement analytically
opt = find_best_focus(res.detector)
println("Nominal detector y: ", opt.y_nom * 1e3, " mm")
println("Optimal detector y: ", opt.y_opt * 1e3, " mm (Δy = ", opt.delta_y * 1e6, " µm)")
println("RMS nominal: ", opt.rms_nom * 1e6, " µm -> RMS optimal: ", opt.rms_opt * 1e6, " µm")

# Shift detector in-place to optimal focus
refocus!(res)
```

### 4. Generate Standalone Julia Code from ZMX

Convert any `.zmx` file into an independent Julia script:

```julia
using BeamletOpticsPrescriptionIO

generate_bmo_script("input_lens.zmx", "standalone_model.jl")
```

The resulting script contains pure, idiomatic `BeamletOptics.jl` constructors:

```julia
using BeamletOptics, LinearAlgebra
const mm = 1e-3

lens_1 = SphericalLens(14.16mm, 69.41mm, 1.386mm, 12.6mm, λ -> 1.5725)
translate3d!(lens_1, [0.0, 4.5mm, 0.0])

doublet_2 = SphericalDoubletLens(16.2mm, -19.66mm, 8.793mm, 2.313mm, 1.386mm, 11.88mm, λ -> 1.5725, λ -> 1.5796)
translate3d!(doublet_2, [0.0, 6.003mm, 0.0])

detector = Detector(45.23mm)
translate3d!(detector, [0.0, 61.88mm, 0.0])

system = System([lens_1, doublet_2, detector])
```

---

## Command-Line Interface (CLI)

A command-line script is provided in `bin/zmx2bmo.jl`:

```bash
# Print system metadata, surfaces, and identified optical elements:
julia --project=. bin/zmx2bmo.jl path/to/lens.zmx --info

# Convert a .zmx file to standalone Julia code:
julia --project=. bin/zmx2bmo.jl path/to/lens.zmx output_model.jl
```

---

## Reference Test Data & Benchmarking

To ensure full reproducibility without committing large third-party proprietary files, a standalone download helper is included.

### Setting Up Reference Data

Run the cross-platform setup script to download external reference databases (with disclaimer prompt):

```bash
julia --project=. scripts/download_test_data.jl
```

This sets up:
1. **Dan Reiley Optical Patent Database**: ~980 optical designs across 9 categories (Endoscopes, Eyepieces, Microscope Objectives, Photographic Primes, Zooms, Projectors, Scan Lenses, Spectrometers, Telescopes).
2. **Thorlabs Zemax Catalog**: 4,217 catalog lenses.

*(All test data directories are automatically ignored by git).*

### Running Validation Benchmarks

```bash
# Run the complete test suite (unit tests & end-to-end raytracing):
julia --project=. test/runtests.jl

# Run the Thorlabs catalog benchmark (evaluates spot size vs. Airy disk):
julia --project=. examples/benchmark_thorlabs.jl --category=achromats

# Run the 9-category patent benchmark:
julia --project=. examples/benchmark_focus_metric.jl --max-per-cat=10
```

---

## Attribution & Notices

Parts of the Zemax token parsing logic and glass name normalization are inspired by the open-source Python packages [`ray-optics`](https://github.com/mjhoptics/ray-optics) and [`opticalglass`](https://github.com/mjhoptics/opticalglass) by **Michael J. Hayford** (licensed under the BSD 3-Clause License). See [NOTICE.md](NOTICE.md) for full license text.

---

## Legal Disclaimer

* **Trademarks**: Zemax® and OpticStudio® are registered trademarks of Zemax, LLC (an Ansys company). Thorlabs® is a registered trademark of Thorlabs, Inc. BeamletOpticsPrescriptionIO.jl is an independent open-source project and is not affiliated with, endorsed by, or sponsored by Zemax, Ansys, or Thorlabs.
* **Third-Party Data**: This repository does not host or redistribute proprietary optical design files. Users are responsible for complying with the terms of service and copyright laws governing any external reference files downloaded.
* For full legal terms, see [DISCLAIMER.md](DISCLAIMER.md).

---

## License

BeamletOpticsPrescriptionIO.jl is released under the [MIT License](LICENSE).
