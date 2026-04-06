# Atomistic Modeling of Laser Powder Bed Fusion Process

<img width="1600" height="1200" alt="lpbf" src="https://github.com/user-attachments/assets/fc29108a-6b28-448c-9158-336d3be096bf" />

[![LAMMPS](https://img.shields.io/badge/LAMMPS-Molecular%20Dynamics-blue)](https://www.lammps.org/)
[![Language](https://img.shields.io/badge/Language-LAMMPS%20Script-orange)](https://docs.lammps.org/)
[![Material](https://img.shields.io/badge/Material-Aluminum%20(Al)-silver)](https://en.wikipedia.org/wiki/Aluminium)
[![GitHub](https://img.shields.io/badge/GitHub-akshansh11-black?logo=github)](https://github.com/akshansh11)

A LAMMPS input script for simulating **Laser Powder Bed Fusion (LPBF)** of Aluminum at the atomistic scale. The script models laser melting of Al powder particles on an Al substrate, including equilibration, laser scanning, and cooling stages — with robust fixes for substrate stability.

---

## Table of Contents

- [Overview](#overview)
- [Root Causes Fixed](#root-causes-fixed)
- [Simulation Stages](#simulation-stages)
- [System Architecture](#system-architecture)
- [Parameters](#parameters)
- [Requirements](#requirements)
- [Usage](#usage)
- [Outputs](#outputs)
- [Visualization](#visualization)

---

## Overview

This simulation models the LPBF process for Aluminum using classical molecular dynamics. A substrate is built from FCC Al, and 15 powder spheres (3×5 grid) are deposited above it. A moving laser then scans across the powder, melting and fusing it to the substrate, followed by a controlled cooling stage.

The script is designed to **prevent substrate atom blowup** — a common failure mode in LPBF MD simulations — through careful thermostat layering and boundary condition management.

---

## Root Causes Fixed

Previous versions of this simulation suffered from substrate atoms being ejected ("blown") due to five identified root causes, all of which are corrected here:

| # | Root Cause | Fix Applied |
|---|-----------|-------------|
| 1 | Langevin damping too slow (`100×dt = 0.1 ps`) — energy built up faster than the thermostat could remove it | Damping reduced to `10×dt = 0.01 ps` (10× faster dissipation) |
| 2 | Only 1 unit cell fixed — insufficient mechanical anchor at the base | Fixed layers increased to **3 unit cells (6 monolayers)** |
| 3 | Only 1 unit cell as thermostat — too thin to absorb conducted heat | Thermostat layers increased to **2 unit cells (4 monolayers)** |
| 4 | `velocity` applied to `g_fixed` group — caused an initial momentum kick to neighbors | Fixed atoms explicitly zeroed: `velocity g_fixed set 0 0 0` |
| 5 | Conduction layers had no energy sink during equilibration | Thermostat velocity rescaled after warm-up to lock in `T_bed` |

---

## Simulation Stages

```
┌─────────────────────────────────────────────────────────────────────┐
│  Stage 0  │  Three-stage Minimization (QuickMin × 2 → CG)          │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 1  │  NVE/limit warm-up (200 steps, displacement cap 0.05 Å) │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 2  │  Equilibration  — 400,000 steps @ 300 K                 │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 3  │  Laser Scan     — steps auto-computed from scan speed    │
├─────────────────────────────────────────────────────────────────────┤
│  Stage 4  │  Cooling        — 2,000,000 steps, NVT ramp to 300 K    │
└─────────────────────────────────────────────────────────────────────┘
```

---

## System Architecture

The substrate is divided into three functional layers stacked in Z:

```
  Z ▲
    │  ┌──────────────────────────┐
    │  │    Powder Spheres (×15)  │  ← 3×5 grid, radius 24.9 Å
    │  │  (NVE, laser target)     │
    │  ├──────────────────────────┤  z_conduct_hi
    │  │   Conduction Layer       │  ← 3 unit cells, NVE
    │  ├──────────────────────────┤  z_therm_hi
    │  │   Thermostat Layer       │  ← 2 unit cells, Langevin @ 300 K
    │  ├──────────────────────────┤  z_fix_hi
    │  │   Fixed Layer            │  ← 3 unit cells, zero force & velocity
    │  └──────────────────────────┘  z = 0
```

**Powder grid:** 3 columns × 5 rows = **15 spheres**, each of radius 24.9 Å, arranged symmetrically above the substrate centre.

---

## Parameters

### Physical

| Parameter | Value | Description |
|-----------|-------|-------------|
| `a0` | 4.05 Å | FCC Al lattice constant |
| `T_bed` | 300 K | Substrate / bed temperature |
| `laser_radius` | 50.0 Å | Gaussian laser spot radius |
| `laser_power` | 0.40 eV/ps | Energy deposition rate |
| `scan_speed` | 0.4 Å/ps | Laser scan velocity (Y direction) |
| `dt` | 0.001 ps | MD timestep |

### Substrate Geometry

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Nx_sub` | 62 UC | X dimension of substrate |
| `Ny_sub` | 111 UC | Y dimension of substrate |
| `Nz_fixed` | 3 UC | Frozen bottom layers |
| `Nz_therm` | 2 UC | Langevin thermostat layers |
| `Nz_conduct` | 3 UC | NVE conduction layers |

### Simulation Length

| Stage | Steps | Notes |
|-------|-------|-------|
| Equilibration | 400,000 | Dump every 500 steps |
| Laser scan | Auto | `ly / scan_speed / dt` |
| Cooling | 2,000,000 | Dump every 100 steps |

---

## Requirements

- **LAMMPS** (any recent stable release, compiled with `MANYBODY` package for EAM)
- **Potential file:** `Al99.eam.alloy` (Mishin EAM potential for Al)
  - Available from the [NIST Interatomic Potentials Repository](https://www.ctcms.nist.gov/potentials/)
- **OVITO** (recommended for visualization)

### LAMMPS Packages Required

```
MANYBODY    # for pair_style eam/alloy
```

---

## Usage

1. **Clone the repository**
   ```bash
   git clone https://github.com/akshansh11/<repo-name>.git
   cd <repo-name>
   ```

2. **Place the potential file** in the same directory as the input script:
   ```
   Al99.eam.alloy
   ```

3. **Run the simulation**
   ```bash
   # Serial
   lammps -in lpbf_al.lammps

   # Parallel (recommended)
   mpirun -np 16 lammps -in lpbf_al.lammps
   ```

4. **Monitor progress** — watch the thermo output columns:
   ```
   Step  Temp  c_T_therm  c_T_pow  PE  KE  Press  Atoms
   ```
   - `c_T_therm` should remain near **300 K** during equilibration
   - `c_T_pow` will spike sharply when the laser passes through

---

## Outputs

| File | Description |
|------|-------------|
| `dump_00_check.lammpstrj` | Initial geometry check (verify spheres in OVITO before full run) |
| `dump_ALL.lammpstrj` | Combined trajectory — all three stages appended |
| `al_pbf_final.data` | Final atomic positions in LAMMPS data format |
| `al_pbf_final.restart` | Restart file for continuing the simulation |

### Per-atom quantities in dump

| Column | Quantity |
|--------|----------|
| `c_cna` | Common Neighbour Analysis (structure identification) |
| `c_pe_atom` | Per-atom potential energy (eV) |
| `c_ke_atom` | Per-atom kinetic energy (eV) |
| `v_T_atom` | Per-atom temperature (K) |
| `c_stress[1–6]` | Per-atom stress tensor components (bar·Å³) |
| `v_S_hyd` | Hydrostatic stress (bar·Å³) |
| `v_S_vm` | Von Mises stress (bar·Å³) |

---

## Visualization

Open trajectory files in **OVITO**:

1. Load `dump_ALL.lammpstrj`
2. Apply **Color Coding** by `v_T_atom` to see the thermal gradient
3. Apply **CNA** modifier — distinguishes FCC (green), HCP (red), and disordered/liquid (white) atoms
4. Use **Slice** modifier to view cross-sections of the melt pool

---

## Author

**akshansh11** — [github.com/akshansh11](https://github.com/akshansh11)

---

## License

This project is open-source. Feel free to use and adapt the script for your own LPBF or MD simulations. Attribution appreciated.
