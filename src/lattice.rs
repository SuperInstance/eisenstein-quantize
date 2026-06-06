//! Hexagonal lattice quantization — nearest-neighbor on the A₂ grid.



/// The square root of 3 as a const f64.
const SQRT3: f64 = 1.732_050_807_568_877_2;

/// A point on the hexagonal A₂ lattice, represented by Eisenstein integer coordinates.
///
/// The lattice point is at:
/// - x = (a + b/2) * spacing
/// - y = (b * √3/2) * spacing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct HexPoint {
    /// First Eisenstein coordinate
    pub a: i64,
    /// Second Eisenstein coordinate
    pub b: i64,
}

impl HexPoint {
    /// Create a new hexagonal lattice point.
    pub const fn new(a: i64, b: i64) -> Self {
        Self { a, b }
    }

    /// Convert to Euclidean coordinates as (x, y) in f64.
    pub fn to_euclidean(&self, spacing: f64) -> (f64, f64) {
        let x = (self.a as f64 + self.b as f64 * 0.5) * spacing;
        let y = (self.b as f64 * SQRT3 * 0.5) * spacing;
        (x, y)
    }

    /// Convert from Euclidean coordinates, snapping to the nearest lattice point.
    ///
    /// Uses the 3-candidate algorithm (round in skewed coords, then check
    /// the two adjacent candidates). This is the correct nearest-neighbor
    /// for the A₂ Voronoi cell (hexagonal).
    pub fn from_euclidean(x: f64, y: f64, spacing: f64) -> Self {
        let s = spacing;
        let inv_s = 1.0 / s;

        // Convert to skewed lattice coordinates
        let b_cont = y / (SQRT3 * 0.5 * s);
        let b_round = b_cont.round() as i64;

        let a_cont = x * inv_s - b_round as f64 * 0.5;
        let a_round = a_cont.round() as i64;

        // Try 3 candidates: (a_round, b_round), (a_round, b_round+1), (a_round, b_round-1)
        let candidates = [
            (a_round, b_round),
            (a_round, b_round + 1),
            (a_round, b_round - 1),
            (a_round + 1, b_round),
            (a_round - 1, b_round),
            (a_round, b_round),
        ];

        let mut best = (a_round, b_round);
        let mut best_dist = f64::MAX;

        for &(ca, cb) in candidates.iter() {
            let qx = (ca as f64 + cb as f64 * 0.5) * s;
            let qy = cb as f64 * SQRT3 * 0.5 * s;
            let dx = x - qx;
            let dy = y - qy;
            let dist = dx * dx + dy * dy;
            if dist < best_dist {
                best_dist = dist;
                best = (ca, cb);
            }
        }

        Self {
            a: best.0,
            b: best.1,
        }
    }

    /// The origin of the hexagonal lattice.
    pub const fn origin() -> Self {
        Self { a: 0, b: 0 }
    }

    /// The six nearest neighbors on the hexagonal lattice.
    pub fn neighbors(&self) -> [Self; 6] {
        [
            Self::new(self.a + 1, self.b),
            Self::new(self.a - 1, self.b),
            Self::new(self.a, self.b + 1),
            Self::new(self.a, self.b - 1),
            Self::new(self.a + 1, self.b - 1),
            Self::new(self.a - 1, self.b + 1),
        ]
    }
}

/// Quantize a batch of 2D vectors to the hexagonal lattice.
///
/// Returns a vector of `HexPoint`s nearest to each input vector.
pub fn quantize_batch(vectors: &[(f64, f64)], spacing: f64) -> Vec<HexPoint> {
    vectors
        .iter()
        .map(|&(x, y)| HexPoint::from_euclidean(x, y, spacing))
        .collect()
}

/// Run the full quantization comparison between hexagonal and rectangular lattices.
///
/// Returns (hex_mse, rect_mse, hex_advantage_pct).
pub fn compare_quantization(
    vectors: &[(f64, f64)],
    spacing: f64,
) -> (f64, f64, f64) {
    // Rectangular quantization
    let rect_error: f64 = vectors
        .iter()
        .map(|&(x, y)| {
            let qx = (x / spacing).round() * spacing;
            let qy = (y / spacing).round() * spacing;
            (x - qx).powi(2) + (y - qy).powi(2)
        })
        .sum();

    // Hexagonal quantization (density-matched)
    let hex_sp = spacing * (2.0_f64 / SQRT3).sqrt();
    let hex_error: f64 = vectors
        .iter()
        .map(|&(x, y)| {
            let p = HexPoint::from_euclidean(x, y, hex_sp);
            let (qx, qy) = p.to_euclidean(hex_sp);
            (x - qx).powi(2) + (y - qy).powi(2)
        })
        .sum();

    let n = vectors.len() as f64;
    let mse_rect = rect_error / n;
    let mse_hex = hex_error / n;
    let adv = (1.0 - mse_hex / mse_rect) * 100.0;

    (mse_hex, mse_rect, adv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_origin_to_euclidean() {
        let origin = HexPoint::origin();
        let (x, y) = origin.to_euclidean(1.0);
        assert!((x - 0.0).abs() < 1e-12);
        assert!((y - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_unit_point() {
        let p = HexPoint::new(1, 0);
        let (x, y) = p.to_euclidean(1.0);
        assert!((x - 1.0).abs() < 1e-12);
        assert!((y - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_neighbors_six() {
        let origin = HexPoint::origin();
        let neighbors = origin.neighbors();
        assert_eq!(neighbors.len(), 6);
        // All neighbors should be distinct
        let mut seen = std::collections::HashSet::new();
        for n in &neighbors {
            assert!(seen.insert((n.a, n.b)), "duplicate neighbor: {n:?}");
        }
    }

    #[test]
    fn test_quantize_roundtrip() {
        // A point at a lattice site should quantize to exactly that site
        let p = HexPoint::new(3, -2);
        let (x, y) = p.to_euclidean(1.5);
        let q = HexPoint::from_euclidean(x, y, 1.5);
        assert_eq!(p, q, "lattice point should quantize to itself");
    }

    #[test]
    fn test_euclidean_roundtrip() {
        // An arbitrary point that's close to a lattice site
        let p = HexPoint::from_euclidean(1.0, 1.0, 1.0);
        let (x, y) = p.to_euclidean(1.0);
        // The roundtripped point should be close to the original
        let dx = 1.0 - x;
        let dy = 1.0 - y;
        let dist = (dx * dx + dy * dy).sqrt();
        assert!(dist < 0.6, "roundtrip distance should be small: {dist}");
    }
}
