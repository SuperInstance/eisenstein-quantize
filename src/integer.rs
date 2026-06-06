//! Eisenstein integer arithmetic — the natural number system for ternary.
//!
//! An Eisenstein integer is `a + bω` where `ω = (-1 + √-3)/2`.
//! These form a Euclidean domain, so they support division with remainder.

use std::ops::{Add, Sub, Mul, Neg};
/// ω = (-1 + √-3)/2 = e^(2πi/3), the primitive cube root of unity.
/// The key property: ω² + ω + 1 = 0, so ω² = -1 - ω.

/// An Eisenstein integer: a + b*ω, where ω = (-1 + √-3)/2.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Eisenstein {
    pub a: i64,
    pub b: i64,
}

impl Eisenstein {
    pub const fn new(a: i64, b: i64) -> Self {
        Self { a, b }
    }

    pub const fn zero() -> Self {
        Self { a: 0, b: 0 }
    }

    pub const fn one() -> Self {
        Self { a: 1, b: 0 }
    }

    /// The unit ω (omega) — primitive cube root of unity.
    pub const fn omega() -> Self {
        Self { a: 0, b: 1 }
    }

    /// Conjugate: a + b*ω̄ = a - b - b*ω
    pub fn conj(&self) -> Self {
        Self {
            a: self.a - self.b,
            b: -self.b,
        }
    }

    /// Norm: N(a + bω) = a² - ab + b²
    pub fn norm(&self) -> i64 {
        self.a * self.a - self.a * self.b + self.b * self.b
    }

    /// Check if this is a unit (norm == 1).
    /// The six units are: ±1, ±ω, ±ω²
    pub fn is_unit(&self) -> bool {
        self.norm() == 1
    }
}

impl Add for Eisenstein {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b,
        }
    }
}

impl Sub for Eisenstein {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            a: self.a - rhs.a,
            b: self.b - rhs.b,
        }
    }
}

impl Neg for Eisenstein {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            a: -self.a,
            b: -self.b,
        }
    }
}

impl Mul for Eisenstein {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // (a₁ + b₁ω)(a₂ + b₂ω) = a₁a₂ + (a₁b₂ + a₂b₁)ω + b₁b₂ω²
        // Since ω² = -1 - ω:
        // = (a₁a₂ - b₁b₂) + (a₁b₂ + a₂b₁ - b₁b₂)ω
        Self {
            a: self.a * rhs.a - self.b * rhs.b,
            b: self.a * rhs.b + self.b * rhs.a - self.b * rhs.b,
        }
    }
}

/// Find all Eisenstein integers within a given norm bound.
pub fn find_by_norm(max_norm: i64) -> Vec<Eisenstein> {
    let bound = (max_norm as f64).sqrt() as i64 + 1;
    let mut result = Vec::new();
    for a in -bound..=bound {
        for b in -bound..=bound {
            let z = Eisenstein::new(a, b);
            if z.norm() <= max_norm && z.norm() > 0 {
                result.push(z);
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_omega_property() {
        // ω satisfies ω² + ω + 1 = 0
        let omega = Eisenstein::omega();
        let omega_sq = omega * omega;
        let sum = omega_sq + omega + Eisenstein::one();
        assert_eq!(sum, Eisenstein::zero(), "ω² + ω + 1 = 0");
    }

    #[test]
    fn test_omega_cubed() {
        let omega = Eisenstein::omega();
        let omega_3 = omega * omega * omega;
        assert_eq!(omega_3, Eisenstein::one(), "ω³ = 1");
    }

    #[test]
    fn test_norm_multiplicative() {
        let z1 = Eisenstein::new(3, 2);
        let z2 = Eisenstein::new(1, -4);
        let prod = z1 * z2;
        assert_eq!(
            z1.norm() * z2.norm(),
            prod.norm(),
            "N(z1 * z2) = N(z1) * N(z2)"
        );
    }

    #[test]
    fn test_six_units() {
        let one = Eisenstein::one();
        let omega = Eisenstein::omega();
        let omega_sq = omega * omega;

        let units = [one, -one, omega, -omega, omega_sq, -omega_sq];
        for u in &units {
            assert!(u.is_unit(), "{u:?} should be a unit");
        }
    }

    #[test]
    fn test_conjugate_property() {
        let z = Eisenstein::new(3, 2);
        let conj = z.conj();
        // z * conj(z) = N(z)
        let product = z * conj;
        assert_eq!(product.a, z.norm());
        assert_eq!(product.b, 0);
    }
}
