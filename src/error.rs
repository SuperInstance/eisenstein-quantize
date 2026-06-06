//! Error types for Eisenstein quantization.

use std::fmt;

/// Errors that can occur during Eisenstein operations.
#[derive(Debug, Clone, PartialEq)]
pub enum EisensteinError {
    /// Vector dimension mismatch (expected 2D)
    DimensionMismatch { expected: usize, got: usize },
    /// Invalid spacing (must be positive)
    InvalidSpacing(f64),
    /// Integer overflow during arithmetic
    Overflow,
}

impl fmt::Display for EisensteinError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch { expected, got } => {
                write!(f, "expected {expected} dimensions, got {got}")
            }
            Self::InvalidSpacing(s) => write!(f, "invalid spacing: {s} (must be positive)"),
            Self::Overflow => write!(f, "integer overflow"),
        }
    }
}

impl std::error::Error for EisensteinError {}
