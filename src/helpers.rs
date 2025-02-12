// helpers.rs
use crate::{AffinePoint, ProjectivePoint};

pub fn show_affine_point(affine_point: &AffinePoint) {
    println!("Affine Point 1: ({}, {})", affine_point.x, affine_point.y)
}

pub fn show_projective_point(projective_point: &ProjectivePoint) {
    println!(
        "Projective Point 1: ({}, {}, {})",
        projective_point.x, projective_point.y, projective_point.z
    )
}
