# Ternary Field — Scalar Field Operations on Ternary Grids

**Ternary Field** computes vector calculus operations — gradient, Laplacian, divergence, and curl — on 2D grids where each cell holds a ternary value {-1, 0, +1}. These differential operators are the foundation for modeling spatial dynamics in ternary systems, from Ising-like lattice models to wave propagation and diffusion.

## Why It Matters

Spatial dynamics matter. A ternary grid where cells hold {-1, 0, +1} values can represent agent territories, temperature states, or material phases. Computing gradients reveals where transitions occur; the Laplacian identifies curvature and diffusion pressure; divergence shows sources and sinks. These are the exact operations needed to model wave propagation, pattern formation, and critical phenomena on ternary lattices. Without differential operators, spatial ternary analysis is limited to counting and neighborhood lookups — with them, you can solve PDEs, detect phase boundaries, and model reaction-diffusion dynamics.

## How It Works

### Gradient

The gradient at cell (x, y) approximates the first derivative:

```
∇f = (∂f/∂x, ∂f/∂y) ≈ (f(x+1,y) - f(x,y), f(x,y+1) - f(x,y))
```

Forward differences at boundaries. Since values are in {-1, 0, +1}, the gradient components are in {-2, -1, 0, 1, 2}. Computed in O(1) per cell, O(N) for the grid.

### Laplacian

The discrete Laplacian measures second-order curvature:

```
Δf = Σ_neighbors f(neighbor) - k · f(center)
```

where k is the number of neighbors (2 at corners, 3 at edges, 4 interior). This is the discrete version of ∂²f/∂x² + ∂²f/∂y². The Laplacian identifies:
- **Positive**: local minimum (concave up)
- **Negative**: local maximum (concave down)
- **Zero**: flat or saddle point

O(N) for the full grid.

### Divergence

Divergence of the gradient field: `∇·(∇f) = ∂²f/∂x² + ∂²f/∂y²`. This equals the Laplacian for a scalar field but is computed differently — as the derivative of each gradient component. Sources (positive divergence) indicate cells where "material" is generated; sinks (negative) indicate absorption.

### Curl

For a 2D scalar field treated as the z-component of a vector potential, curl measures rotation:

```
curl = ∂Fy/∂x - ∂Fx/∂y
```

This identifies rotational flow patterns in the ternary field.

## Quick Start

```rust
use ternary_field::{gradient, laplacian, divergence};

let grid: Vec<i8> = vec![
    0,  0,  0,  0,
    0,  1,  1,  0,
    0,  1, -1,  0,
    0,  0,  0,  0,
];
let width = 4;

let grad = gradient(&grid, width);
let lap = laplacian(&grid, width);
let div = divergence(&grid, width);

// Laplacian at center cells reveals the phase boundary
```

```bash
cargo add ternary-field
```

## API

| Type / Function | Description |
|---|---|
| `gradient(grid, width) → Vec<(f64, f64)>` | ∂/∂x, ∂/∂y per cell |
| `laplacian(grid, width) → Vec<f64>` | Second derivative (curvature) |
| `divergence(grid, width) → Vec<f64>` | Sources and sinks |

## Architecture Notes

Field operators enable spatial reasoning in **SuperInstance**: detecting phase boundaries, modeling agent territory dynamics, and computing pressure gradients across fleet nodes. The γ + η = C conservation law has a field-theoretic interpretation: the divergence of the γ field (growth flux) must equal the η sink rate (entropy absorption) at every point. See [Architecture](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

- Strang, Gilbert. *Computational Science and Engineering*, Wellesley-Cambridge, 2007 — finite difference methods.
- Press, W. H. et al. *Numerical Recipes*, 3rd ed., Cambridge UP, 2007 — grid-based PDE solvers.
- Trefethen, Lloyd N. *Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations*, Cornell, 1996.

## License

MIT
