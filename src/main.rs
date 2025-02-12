mod helpers;

use num_bigint::BigInt;
use num_traits::{Num, One, Pow, Zero};

#[derive(Clone, Debug)]
struct AffinePoint {
    x: BigInt,
    y: BigInt,
}

#[derive(Clone, Debug)]
struct ProjectivePoint {
    x: BigInt,
    y: BigInt,
    z: BigInt,
}

struct Curve {
    p: BigInt,
    a: BigInt,
    b: BigInt,
}

impl AffinePoint {
    fn new(x: BigInt, y: BigInt) -> AffinePoint {
        AffinePoint { x, y }
    }

    fn to_projective(&self) -> ProjectivePoint {
        ProjectivePoint {
            x: self.x.clone(),
            y: self.y.clone(),
            z: BigInt::zero(),
        }
    }
}

impl ProjectivePoint {
    fn new(x: BigInt, y: BigInt, z: BigInt) -> Self {
        ProjectivePoint { x, y, z }
    }

    fn neutral() -> Self {
        ProjectivePoint {
            x: BigInt::zero(),
            y: BigInt::one(),
            z: BigInt::zero(),
        }
    }

    fn to_affine(&self, curve: &Curve) -> Option<AffinePoint> {
        match self.z.clone().modinv(&curve.p) {
            Some(z_inv) => {
                let x: BigInt = (&self.x * &z_inv) % &curve.p;
                let y: BigInt = (&self.y * &z_inv) % &curve.p;

                Some(AffinePoint::new(x, y))
            }
            None => None,
        }
    }
}

impl Curve {
    fn new(p: &BigInt, a: &BigInt, b: &BigInt) -> Option<Curve> {
        let d: BigInt = &a.pow(3u32) * 4 + &b.pow(2u32) * 27;

        match d != BigInt::zero() {
            true => Some(Curve {
                p: p.clone(),
                a: a.clone(),
                b: b.clone(),
            }),
            false => None,
        }
    }

    fn verify_affine_point(&self, point: &AffinePoint) -> bool {
        let left: BigInt = point.clone().y.pow(2u32);
        let right: BigInt = point.clone().x.pow(3u32) + &self.a * point.clone().x + &self.b;

        left == right
    }

    fn verify_projective_point(&self, point: &ProjectivePoint) -> bool {
        let z1: BigInt = point.z.clone();
        let z2: BigInt = point.z.clone().pow(2u32);
        let z3: BigInt = point.z.clone().pow(3u32);

        let left: BigInt = point.clone().y.pow(2u32) * z1;
        let right: BigInt =
            point.clone().x.pow(3u32) + &self.a * z2 * point.clone().x + &self.b * z3;

        left == right
    }

    fn add(&self, point1: &ProjectivePoint, point2: &ProjectivePoint) -> ProjectivePoint {
        ProjectivePoint::neutral()
    }

    fn double(&self) -> ProjectivePoint {
        ProjectivePoint::neutral()
    }

    fn scalar_mul(&self, scalar: &BigInt) -> ProjectivePoint {
        ProjectivePoint::neutral()
    }
}

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_to_affine_point_trivial() {
        let p: &BigInt = &BigInt::from_str_radix("26", 10).unwrap();
        let a: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();
        let b: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 16).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 16).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), BigInt::one());

        // (2, 3)
        let affine_point: AffinePoint = projective_point.to_affine(&curve).unwrap();

        assert_eq!(&affine_point.x, &x);
        assert_eq!(&affine_point.y, &y);
    }

    #[test]
    fn test_map_to_affine_point_inverse_exist() {
        let p: &BigInt = &BigInt::from_str_radix("26", 10).unwrap();
        let a: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();
        let b: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 10).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 10).unwrap();
        let z: BigInt = BigInt::from_str_radix("11", 10).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        // (12, 5)
        let affine_point: AffinePoint = projective_point.to_affine(&curve).unwrap();

        assert_eq!(affine_point.x, BigInt::from(12));
        assert_eq!(affine_point.y, BigInt::from(5));
    }

    #[test]
    fn test_map_to_affine_point_inverse_non_exist() {
        let p: &BigInt = &BigInt::from_str_radix("26", 10).unwrap();
        let a: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();
        let b: &BigInt = &BigInt::from_str_radix("1", 10).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("2", 10).unwrap();
        let y: BigInt = BigInt::from_str_radix("3", 10).unwrap();
        let z: BigInt = BigInt::from_str_radix("10", 10).unwrap();

        // (2, 3, 1)
        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        // (12, 5)
        let affine_point: Option<AffinePoint> = projective_point.to_affine(&curve);
        assert!(&affine_point.is_none());
    }

    // P-256
    #[test]
    fn test_add_affine_point() {
        let p: &BigInt = &BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
            16,
        )
        .unwrap();
        let a: &BigInt = &BigInt::from_str_radix(
            "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
            16,
        )
        .unwrap();
        let b: &BigInt = &BigInt::from_str_radix(
            "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            16,
        )
        .unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();

        let x1: BigInt = BigInt::from_str_radix("2", 16).unwrap();
        let y1: BigInt = BigInt::from_str_radix("3", 16).unwrap();
        let z1: BigInt = BigInt::from_str_radix("4", 16).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let identity: ProjectivePoint = ProjectivePoint::neutral();

        let result: ProjectivePoint = curve.add(&point1, &identity);
    }
}
