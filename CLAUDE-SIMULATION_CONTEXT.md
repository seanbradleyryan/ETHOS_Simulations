# SIMULATION_CONTEXT.md — k-Wave, Tissue Models & Sensor Design

> Referenced when working on: `run_single_field_simulation.m`, `run_standalone_simulation.m`, `create_acoustic_medium.m`, `determine_sensor_mask.m`, `apply_element_averaging.m`, `find_optimal_kwave_size.m`, `run_medium_comparison.m`, `test_time_dependence.m`

## k-Wave Simulation Parameters

| Parameter | Value |
|-----------|-------|
| PML size | 10 voxels (default) |
| CFL number | 0.3 |
| Time reversal BCs | Dirichlet boundary conditions |
| GPU acceleration | `use_gpu = true` |
| Typical runtime | ~23 min for 7 fields on 256×256×128 grid with GPU |

## Tissue / Grüneisen Methods (`sensor_placement_method`)

| Method | Description |
|--------|-------------|
| `'uniform'` | Single properties everywhere: Γ=1.0, c=1540 m/s, ρ=1000 kg/m³ |
| `'threshold_1'` | 9 tissue types by HU: air, lung, fat, water, blood, muscle, soft tissue, bone, metal |
| `'threshold_2'` | 4 tissue types by HU: water, fat, soft tissue, bone — **most commonly used** |

## Sensor Design

- **Geometry:** Flat rigid 10×10 cm planar array (not curved). Real arrays are rigid; tissue deforms, not the sensor.
- **Placement:** Anterior abdomen, avoiding all beam field jaw projections on the anterior surface.
- **Coupling:** Water fills gaps between sensor and body; no coupling concerns outside body.
- **Beam exclusion:** Uses divergence geometry to project jaw openings from isocenter onto anterior surface.

### Sensor Placement Methods (`sensor_placement_method`)

| Value | Description |
|-------|-------------|
| `'full_plane_anterior'` | Full YZ plane at `sensor_x_index` *(default)* |
| `'full_plane_lateral'` | Full XZ plane at `sensor_y_index` |
| `'spherical'` | Spherical shell geometry |

## Gotchas

- **Slice alignment streaking:** Misaligned CT/dose grids cause horizontal streaking artifacts in body masks. Always verify grid alignment after resampling.
- **GPU memory:** 256³ grids with complex tissue models can exhaust GPU memory. Monitor with `gpuDevice`; fall back to CPU if needed.
