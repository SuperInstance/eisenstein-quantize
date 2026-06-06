//! # Eisenstein Quantize — Hexagonal A₂ Lattice
//!
//! The densest 2D lattice packing (Thue's theorem, 1890). The hexagonal grid
//! is the natural coordinate system for balanced ternary — Eisenstein integers
//! ℤ[ω] form the A₂ root lattice, which is the 2D projection of the E₈ lattice
//! family used in topological quantum computing.
//!
//! ## Why hexagons are better
//!
//! | Metric | Square Z² | Hexagonal A₂ | Advantage |
//! |--------|-----------|--------------|-----------|
//! | Packing density | π/4 ≈ 0.785 | π/(2√3) ≈ 0.907 | **+15.5%** |
//! | Quantization MSE | baseline | ~3.9% lower | **-3.9%** |
//! | Nearest neighbors | 4 | 6 | **+50%** |
//! | Symmetry group | D₄ (order 8) | D₆ (order 12) | **+50%** |
//!
//! ## Eisenstein integers
//!
//! An Eisenstein integer is `a + bω` where `ω = (-1 + √-3)/2` is the primitive
//! cube root of unity. This maps directly to hexagonal coordinates (a, b) where:
//! - x = (a + b/2) * spacing
//! - y = (b * √3/2) * spacing

pub mod lattice;
pub mod integer;
pub mod error;
