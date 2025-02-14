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
        if point1.x == BigInt::zero() && point1.y == BigInt::one() && point1.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        if point2.x == BigInt::zero() && point2.y == BigInt::one() && point2.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        // U1 := Y2*Z1
        let u1: BigInt = point2.y.clone() * point1.z.clone();

        // U2 := Y1*Z2
        let u2: BigInt = point1.y.clone() * point2.z.clone();

        // V1 := X2*Z1
        let v1: BigInt = point2.x.clone() * point1.z.clone();

        // V2 := X1*Z2
        let v2: BigInt = point1.x.clone() * point2.z.clone();

        if v1 == v2 {
            return if u1 != u2 {
                ProjectivePoint::neutral()
            } else {
                self.double(&point1)
            };
        }

        // U := U1 - U2
        let u: BigInt = u1.clone() - u2.clone();

        // V := V1 - V2
        let v: BigInt = v1.clone() - v2.clone();

        // W := Z1*Z2
        let w: BigInt = point1.z.clone() * point2.z.clone();

        // A := U^2*W - V^3 - 2*V^2*V2
        let a: BigInt = u.clone() * u.clone() * w.clone()
            - v.clone() * v.clone() * v.clone()
            - 2 * v.clone() * v.clone() * v2.clone();

        // X3 := V*A
        let x: BigInt = v.clone() * a.clone();

        // Y3 := U*(V^2*V2 - A) - V^3*U2
        let y: BigInt = u.clone() * (v.clone() * v.clone() * v2.clone() - a.clone())
            - v.clone() * v.clone() * v.clone() * u2.clone();

        // Z3 := V^3*W
        let z: BigInt = v.clone() * v.clone() * v.clone() * w.clone();

        ProjectivePoint::new(x % 29, y % 29, z % 29)
    }

    fn double(&self, point: &ProjectivePoint) -> ProjectivePoint {
        if point.x == BigInt::zero() && point.y == BigInt::one() && point.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        if point.y == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        let w: BigInt = &self.a * point.z.clone();
        let s: BigInt = point.y.clone() * point.z.clone();
        let b: BigInt = point.x.clone() * point.y.clone() * s.clone();
        let h: BigInt = w.clone() * w.clone() - 8 * b.clone();

        let x: BigInt = 2 * h.clone() * s.clone();
        let y: BigInt = w.clone() * (4 * b.clone() - h.clone())
            - 8 * point.y.clone() * point.y.clone() * s.clone() * s.clone();
        let z: BigInt = 8 * s.clone() * s.clone() * s.clone();

        ProjectivePoint::new(x % 29, y % 29, z % 29)
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

    #[test]
    fn test_double1() {
        let p: &BigInt = &BigInt::from_str_radix("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16).unwrap();
        let a: &BigInt = &BigInt::from_str_radix("ffffffff00000001000000000000000000000000fffffffffffffffffffffffc", 16).unwrap();
        let b: &BigInt = &BigInt::from_str_radix("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16).unwrap();

        let curve: Curve = Curve::new(p, a, b).unwrap();
        let x: BigInt = BigInt::from_str_radix("96f43acb7e4ad862d547681483abfe1c6a21eace101ceca1e0dfa97beb8d99ec", 16).unwrap();
        let y: BigInt = BigInt::from_str_radix("69e14d9fe6ce5b8cb04b2fb0ffa525526c2b7114d891b4aa2f53dfe8de3d6dad", 16).unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 16).unwrap();

        let x2: BigInt = BigInt::from_str_radix("23256713097452873534819684224181198488753197392778987539588939509885686328462", 10).unwrap();
        let y2: BigInt = BigInt::from_str_radix("23560293084035690959730279798706588908809082944968261336868665854561491207411", 10).unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let projective_point2: ProjectivePoint =
            ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        println!("{:?}", curve.double(&projective_point));
        println!(
            "{:?}",
            curve.double(&projective_point).to_affine(&curve).unwrap()
        );

        println!();

        println!("{:?}", curve.add(&projective_point, &projective_point2));
        println!(
            "{:?}",
            curve
                .add(&projective_point, &projective_point2)
                .to_affine(&curve)
                .unwrap()
        );
    }

    // P-256
    #[test]
    fn test_add_affine_point_neutral() {
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

    // #[test]
    // fn test_add_two_affine_points() {
    //     let p: &BigInt = &BigInt::from_str_radix(
    //         "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
    //         16,
    //     )
    //         .unwrap();
    //     let a: &BigInt = &BigInt::from_str_radix(
    //         "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
    //         16,
    //     )
    //         .unwrap();
    //     let b: &BigInt = &BigInt::from_str_radix(
    //         "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
    //         16,
    //     )
    //         .unwrap();
    //
    //     let curve: Curve = Curve::new(p, a, b).unwrap();
    //
    //     let x1: BigInt = BigInt::from_str_radix("96f43acb7e4ad862d547681483abfe1c6a21eace101ceca1e0dfa97beb8d99ec", 16).unwrap();
    //     let y1: BigInt = BigInt::from_str_radix("69e14d9fe6ce5b8cb04b2fb0ffa525526c2b7114d891b4aa2f53dfe8de3d6dad", 16).unwrap();
    //     let z1: BigInt = BigInt::from_str_radix("1", 16).unwrap();
    //
    //     let x2: BigInt = BigInt::from_str_radix("4faeebce6274ce7715c658a430621627fef6ef248b1655121a28550ab0099379", 16).unwrap();
    //     let y2: BigInt = BigInt::from_str_radix("e6c1da595c7e174ee06f298d50261ed7256fdb79e42ebf519137309d042ee0ae", 16).unwrap();
    //     let z2: BigInt = BigInt::from_str_radix("1", 16).unwrap();
    //
    //     let xs: BigInt = BigInt::from_str_radix("60925061547489862940527908037003633368921254513905133864735524490867813970614", 10).unwrap();
    //     let ys: BigInt = BigInt::from_str_radix("4224995115825146765508233942604852388022531386580408859355570965592456321292", 10).unwrap();
    //     let zs: BigInt = BigInt::from_str_radix("1", 10).unwrap();
    //
    //     let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
    //     let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());
    //     let sum_point: ProjectivePoint = ProjectivePoint::new(xs.clone(), ys.clone(), zs.clone());
    //
    //
    //     let result: ProjectivePoint = curve.add(&point1, &point2);
    //
    //
    //     println!("{:?}", curve.add(&projective_point, &projective_point2));
    //     println!(
    //         "{:?}",
    //         curve
    //             .add(&projective_point, &projective_point2)
    //             .to_affine(&curve)
    //             .unwrap()
    //     );
    // }
}
