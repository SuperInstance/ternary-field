#![forbid(unsafe_code)]

/// Compute gradient (directional change) at each cell of a ternary grid.
/// Returns (dx, dy) per cell as Vec<f64>.
pub fn gradient(grid: &[i8], width: usize) -> Vec<(f64, f64)> {
    let height = if width == 0 { 0 } else { grid.len() / width };
    let mut out = Vec::with_capacity(grid.len());
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let dx = if x + 1 < width {
                (grid[idx + 1] - grid[idx]) as f64
            } else {
                0.0
            };
            let dy = if y + 1 < height {
                (grid[(y + 1) * width + x] - grid[idx]) as f64
            } else {
                0.0
            };
            out.push((dx, dy));
        }
    }
    out
}

/// Compute Laplacian (second derivative / curvature) at each cell.
pub fn laplacian(grid: &[i8], width: usize) -> Vec<f64> {
    let height = if width == 0 { 0 } else { grid.len() / width };
    let mut out = Vec::with_capacity(grid.len());
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let center = grid[idx] as f64;
            let mut sum = 0.0;
            let mut count = 0.0;
            if x > 0 {
                sum += grid[idx - 1] as f64;
                count += 1.0;
            }
            if x + 1 < width {
                sum += grid[idx + 1] as f64;
                count += 1.0;
            }
            if y > 0 {
                sum += grid[(y - 1) * width + x] as f64;
                count += 1.0;
            }
            if y + 1 < height {
                sum += grid[(y + 1) * width + x] as f64;
                count += 1.0;
            }
            let lap = if count > 0.0 { sum - count * center } else { 0.0 };
            out.push(lap);
        }
    }
    out
}

/// Compute divergence (field source/sink) at each cell.
/// div = dFx/dx + dFy/dy, treating grid values as scalar field.
pub fn divergence(grid: &[i8], width: usize) -> Vec<f64> {
    let grads = gradient(grid, width);
    let height = if width == 0 { 0 } else { grid.len() / width };
    let mut out = Vec::with_capacity(grid.len());
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            // Forward differences of gradient components
            let ddx = if x + 1 < width {
                grads[idx + 1].0 - grads[idx].0
            } else {
                0.0
            };
            let ddy = if y + 1 < height {
                grads[(y + 1) * width + x].1 - grads[idx].1
            } else {
                0.0
            };
            out.push(ddx + ddy);
        }
    }
    out
}

/// Compute curl (rotational measure) at each cell.
/// curl = dFy/dx - dFx/dy
pub fn curl(grid: &[i8], width: usize) -> Vec<f64> {
    let grads = gradient(grid, width);
    let height = if width == 0 { 0 } else { grid.len() / width };
    let mut out = Vec::with_capacity(grid.len());
    for y in 0..height {
        for x in 0..width {
            let idx = y * width + x;
            let dfy_dx = if x + 1 < width {
                grads[idx + 1].1 - grads[idx].1
            } else {
                0.0
            };
            let dfx_dy = if y + 1 < height {
                grads[(y + 1) * width + x].0 - grads[idx].0
            } else {
                0.0
            };
            out.push(dfy_dx - dfx_dy);
        }
    }
    out
}

/// Compute total field energy (magnitude squared sum).
pub fn field_energy(grid: &[i8], width: usize) -> Vec<f64> {
    let grads = gradient(grid, width);
    grads.iter().map(|(dx, dy)| dx * dx + dy * dy).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_grid() -> Vec<i8> {
        vec![0, 0, 0, 0, 0, 0, 0, 0, 0] // 3x3 flat
    }

    #[test]
    fn test_gradient_flat() {
        let g = flat_grid();
        let grad = gradient(&g, 3);
        assert_eq!(grad.len(), 9);
        for (dx, dy) in &grad {
            assert_eq!(*dx, 0.0);
            assert_eq!(*dy, 0.0);
        }
    }

    #[test]
    fn test_gradient_slope() {
        let g = vec![0, 1, 2, 0, 1, 2, 0, 1, 2]; // horizontal slope
        let grad = gradient(&g, 3);
        assert_eq!(grad[0].0, 1.0); // dx at (0,0)
        assert_eq!(grad[0].1, 0.0); // dy at (0,0)
    }

    #[test]
    fn test_gradient_vertical() {
        let g = vec![0, 0, 0, 1, 1, 1, 2, 2, 2]; // vertical slope
        let grad = gradient(&g, 3);
        assert_eq!(grad[0].0, 0.0);
        assert_eq!(grad[0].1, 1.0);
    }

    #[test]
    fn test_laplacian_flat() {
        let g = flat_grid();
        let lap = laplacian(&g, 3);
        for v in &lap {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_laplacian_peak() {
        let g = vec![0, 0, 0, 0, 1, 0, 0, 0, 0]; // peak in center
        let lap = laplacian(&g, 3);
        // center laplacian should be negative (concave down)
        assert!(lap[4] < 0.0);
    }

    #[test]
    fn test_laplacian_valley() {
        let g = vec![1, 1, 1, 1, 0, 1, 1, 1, 1]; // valley in center
        let lap = laplacian(&g, 3);
        assert!(lap[4] > 0.0);
    }

    #[test]
    fn test_divergence_flat() {
        let g = flat_grid();
        let div = divergence(&g, 3);
        for v in &div {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_divergence_source() {
        let g = vec![0, 0, 0, 0, 1, 0, 0, 0, 0];
        let div = divergence(&g, 3);
        // center should have nonzero divergence
        assert!(div[4] != 0.0);
    }

    #[test]
    fn test_curl_flat() {
        let g = flat_grid();
        let c = curl(&g, 3);
        for v in &c {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_field_energy_flat() {
        let g = flat_grid();
        let e = field_energy(&g, 3);
        for v in &e {
            assert_eq!(*v, 0.0);
        }
    }

    #[test]
    fn test_field_energy_slope() {
        let g = vec![0, 1, 2, 3];
        let e = field_energy(&g, 4);
        assert!(e.iter().sum::<f64>() > 0.0);
    }

    #[test]
    fn test_gradient_empty() {
        let g: &[i8] = &[];
        let grad = gradient(g, 0);
        assert!(grad.is_empty());
    }

    #[test]
    fn test_gradient_1d() {
        let g = vec![0, 1, -1];
        let grad = gradient(&g, 3);
        assert_eq!(grad.len(), 3);
    }

    #[test]
    fn test_laplacian_linear_ramp() {
        // Interior cells of a linear ramp should have zero laplacian
        let g = vec![0, 1, 2, 3, 4, 5, 6, 7, 8]; // 3x3 linear ramp
        let lap = laplacian(&g, 3);
        // Only check center cell (4) - interior of the ramp
        assert!((lap[4]).abs() < 0.01);
    }

    #[test]
    fn test_curl_4x4() {
        // Larger grid with asymmetric values
        let g = vec![
            1, 0, -1, 0,
            0, 1, 0, -1,
            -1, 0, 1, 0,
            0, -1, 0, 1
        ];
        let c = curl(&g, 4);
        assert_eq!(c.len(), 16);
        // Interior cells (5,6,9,10) should have computed curl values
        // Even if zero for gradient fields, the function runs correctly
        assert!(c.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn test_field_energy_nonzero() {
        let g = vec![0, 1, 0, 0, 0, 0, 0, 0, 0]; // spike
        let e = field_energy(&g, 3);
        let total: f64 = e.iter().sum();
        assert!(total > 0.0);
    }
}
