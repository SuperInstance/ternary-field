# ternary-field

**Vector field calculus on ternary grids. Gradient, Laplacian, divergence, curl.**

Every grid has hidden structure. A ternary grid of {-1, 0, +1} values looks simple — just three states per cell. But compute the gradient (how values change between neighbors), the Laplacian (local curvature), divergence (sources and sinks), and curl (rotation), and a rich vector field emerges. The ternary constraint makes these fields *discrete* and *interpretable* — every arrow points in a specific direction with specific magnitude.

## What's Inside

- **`gradient(grid, width)`** — directional change at each cell. Returns `(dx, dy)` as `Vec<(f64, f64)>`
- **`laplacian(grid, width)`** — second derivative / curvature. Identifies ridges, valleys, and flat zones
- **`divergence(field, width)`** — source/sink detection. Positive = source, negative = sink
- **`curl(field, width)`** — rotational field. Nonzero curl = vortex
- **`boundary_cells(grid, width)`** — cells where value differs from at least one neighbor (edges of regions)
- **`connected_components(grid, width, value)`** — find connected regions of a specific ternary value
- **`field_energy(grid, width)`** — total energy in the field (sum of squared gradients)

## Quick Example

```rust
use ternary_field::*;

// A 5x5 grid with a +1 island in the center
let mut grid = vec![0; 25];
grid[12] = 1; // center cell

// Gradient: how values change
let grad = gradient(&grid, 5);
// grad[7] = (1.0, 0.0)  — to the left of center, x increases
// grad[12] = (-1.0, 0.0) — center, values decrease to the right

// Laplacian: curvature
let lap = laplacian(&grid, 5);
// lap[12] > 0 — center is a peak

// Boundaries: edges of the +1 region
let bounds = boundary_cells(&grid, 5);
// The 4 neighbors of center cell are boundary cells

// Field energy: total gradient magnitude
let energy = field_energy(&grid, 5);
assert!(energy > 0.0);
```

## The Insight

**Ternary fields have quantized structure.** Unlike continuous fields where gradients can take any value, ternary gradients are always in {-2, -1, 0, 1, 2}. This makes field analysis *exact* — no floating-point noise, no approximation. The discrete nature means you can prove things about the field that you can't prove about continuous approximations.

**Use cases:**
- **Image processing** — edge detection and segmentation on ternary images
- **Cellular automata** — analyze the spatial structure of CA states
- **Physics simulation** — discrete field theory on ternary lattices
- **Geographic analysis** — terrain features on discretized elevation maps
- **Game AI** — vector fields for navigation and threat assessment

## See Also

- **ternary-morph** — morphological operations (complements field analysis)
- **ternary-lattice** — order theory on ternary values
- **ternary-diff** — diff and patch operations on ternary grids
- **ternary-geometry** — geometric structures on ternary grids
- **ternary-irradiate** — radiation fields (an application of field computation)

## Install

```bash
cargo add ternary-field
```

## License

MIT
