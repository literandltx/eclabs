mod helpers;

use crate::helpers::module;
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
            z: BigInt::one(),
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

    fn add(&self, a: &ProjectivePoint, b: &ProjectivePoint) -> ProjectivePoint {
        if a.x == BigInt::zero() && a.y == BigInt::one() && a.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        if b.x == BigInt::zero() && b.y == BigInt::one() && b.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        // U1 := Y2*Z1
        let u1: BigInt = (&b.y * &a.z) % &self.p;

        // U2 := Y1*Z2
        let u2: BigInt = (&a.y * &b.z) % &self.p;

        // V1 := X2*Z1
        let v1: BigInt = (&b.x * &a.z) % &self.p;

        // V2 := X1*Z2
        let v2: BigInt = (&a.x * &b.z) % &self.p;

        if v1 == v2 {
            return if u1 != u2 {
                ProjectivePoint::neutral()
            } else {
                println!("called double");
                self.double(&a)
            };
        }

        // U := U1 - U2
        let u: BigInt = (&u1 - &u2) % &self.p;

        // V := V1 - V2
        let v: BigInt = (&v1 - &v2) % &self.p;

        // W := Z1*Z2
        let w: BigInt = (&a.z * &b.z) % &self.p;

        // A := U^2*W - V^3 - 2*V^2*V2
        let _a: BigInt = (&u * &u * &w - &v * &v * &v - 2 * &v * &v * &v2) % &self.p;

        // X3 := V*A
        let x: BigInt = (&v * &_a) % &self.p;

        // Y3 := U*(V^2*V2 - A) - V^3*U2
        let y: BigInt = (&u * (&v * &v * &v2 - &_a) - &v * &v * &v * &u2) % &self.p;

        // Z3 := V^3*W
        let z: BigInt = (&v * &v * &v * &w) % &self.p;

        ProjectivePoint::new(
            module(&x, &self.p),
            module(&y, &self.p),
            module(&z, &self.p),
        )
    }

    fn double(&self, point: &ProjectivePoint) -> ProjectivePoint {
        if point.x == BigInt::zero() && point.y == BigInt::one() && point.z == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        if point.y == BigInt::zero() {
            return ProjectivePoint::neutral();
        }

        // W := a*Z^2 + 3*X^2
        let w: BigInt = ((&self.a * &point.z * &point.z) + (3 * &point.x * &point.x)) % &self.p;

        // S := Y*Z
        let s: BigInt = (&point.y * &point.z) % &self.p;

        // B := X*Y*S
        let b: BigInt = (&point.x * &point.y * &s) % &self.p;

        // H := W^2 - 8*B
        let h: BigInt = (&w * &w - 8 * &b) % &self.p;

        // X' := 2*H*S
        let x: BigInt = (2 * &h * &s) % &self.p;

        // Y' := W*(4*B - H) - 8*Y^2*S^2
        let y: BigInt = (&w * (4 * &b - &h) - 8 * &point.y * &point.y * &s * &s) % &self.p;

        // Z' := 8*S^3
        let z: BigInt = (8 * &s * &s * &s) % &self.p;

        ProjectivePoint::new(
            module(&x, &self.p),
            module(&y, &self.p),
            module(&z, &self.p),
        )
    }

    fn scalar_mul(&self, scalar: &BigInt, point: &ProjectivePoint) -> ProjectivePoint {
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
        let x: BigInt = BigInt::from_str_radix(
            "96f43acb7e4ad862d547681483abfe1c6a21eace101ceca1e0dfa97beb8d99ec",
            16,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "69e14d9fe6ce5b8cb04b2fb0ffa525526c2b7114d891b4aa2f53dfe8de3d6dad",
            16,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 16).unwrap();

        let double_x: BigInt = BigInt::from_str_radix(
            "23256713097452873534819684224181198488753197392778987539588939509885686328462",
            10,
        )
        .unwrap();
        let double_y: BigInt = BigInt::from_str_radix(
            "23560293084035690959730279798706588908809082944968261336868665854561491207411",
            10,
        )
        .unwrap();
        let double_z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let projective_point: ProjectivePoint =
            ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let double_projective_point: ProjectivePoint = curve.double(&projective_point);
        let double_affine_point = double_projective_point.to_affine(&curve).unwrap();

        let result: ProjectivePoint = double_affine_point.to_projective();

        assert_eq!(result.z, double_z);
        assert_eq!(result.x, double_x);
        assert_eq!(result.y, double_y);
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

    #[test]
    fn test_add_two_affine_points() {
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

        let x1: BigInt = BigInt::from_str_radix(
            "68278443758753309620248260190700068206105694628838492838163118930873086220780",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "47890925436225047886667813019233725009834729938480004609112948108082127138221",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x2: BigInt = BigInt::from_str_radix(
            "36041773901858648762789006722537555576806908678254251775422211831047989597049",
            10,
        )
        .unwrap();
        let y2: BigInt = BigInt::from_str_radix(
            "104374463647532952902541679499731110179928673815908965998191930674877316391086",
            10,
        )
        .unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 16).unwrap();

        let xs: BigInt = BigInt::from_str_radix(
            "60925061547489862940527908037003633368921254513905133864735524490867813970614",
            10,
        )
        .unwrap();
        let ys: BigInt = BigInt::from_str_radix(
            "4224995115825146765508233942604852388022531386580408859355570965592456321292",
            10,
        )
        .unwrap();
        let zs: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        let add_projective_point: ProjectivePoint = curve.add(&point1, &point2);

        let add_affine_point: AffinePoint = add_projective_point.to_affine(&curve).unwrap();

        let resut: ProjectivePoint = add_affine_point.to_projective();

        assert_eq!(resut.x, xs);
        assert_eq!(resut.y, ys);
        assert_eq!(resut.z, zs);
    }

    #[test]
    fn test_multiply_point_on_scalar() {
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

        let x: BigInt = BigInt::from_str_radix(
            "36468809527719095167323546128622644075555946232974359005635481469739685065250",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "55086920727300960012414320103293286937832915861215179864890281344876589042713",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar: BigInt = BigInt::from_str_radix(
            "50313277712876843208627822877019546215941224990452903113298304423968512062362",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "105471300999999635061948442487927318986215238139837487327895353516078990089204",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "102129873108155434225483148806386549747776446537842754889946033852013914809773",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let mul_projective_point: ProjectivePoint = curve.scalar_mul(&scalar, &point);

        let mul_affine_point: AffinePoint = mul_projective_point.to_affine(&curve).unwrap();

        let resut: ProjectivePoint = mul_affine_point.to_projective();

        assert_eq!(resut.x, x_res);
        assert_eq!(resut.y, y_res);
        assert_eq!(resut.z, z_res);
    }

    #[test]
    fn test_scalar_multiplication_using_order() {
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

        let x: BigInt = BigInt::from_str_radix(
            "88830848574649868340568578826653665181106303545668410351759331370193752307474",
            10,
        )
        .unwrap();
        let y: BigInt = BigInt::from_str_radix(
            "22058354116588620035967384506492701987686288267863232015425921192621894353676",
            10,
        )
        .unwrap();
        let z: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar: BigInt = BigInt::from_str_radix(
            "115792089210356248762697446949407573529996955224135760342422259061068512044369",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();
        let y_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();
        let z_res: BigInt = BigInt::from_str_radix("0", 10).unwrap();

        let point: ProjectivePoint = ProjectivePoint::new(x.clone(), y.clone(), z.clone());

        let mul_projective_point: ProjectivePoint = curve.scalar_mul(&scalar, &point);

        let mul_affine_point: AffinePoint = mul_projective_point.to_affine(&curve).unwrap();

        let resut: ProjectivePoint = mul_affine_point.to_projective();

        assert_eq!(resut.x, x_res);
        assert_eq!(resut.y, y_res);
        assert_eq!(resut.z, z_res);
    }

    #[test]
    fn test_complex_scalar_multiplication_and_sum() {
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

        let x1: BigInt = BigInt::from_str_radix(
            "42962812987195471550265028018049668228525471650158073055566585942623576308426",
            10,
        )
        .unwrap();
        let y1: BigInt = BigInt::from_str_radix(
            "85736017906299720575115626781728330074318681691122158879669393468418823605282",
            10,
        )
        .unwrap();
        let z1: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let x2: BigInt = BigInt::from_str_radix(
            "115453654416017434409541019576558067958554154123817093231752062644606305094870115453654416017434409541019576558067958554154123817093231752062644606305094870",
            10,
        )
            .unwrap();
        let y2: BigInt = BigInt::from_str_radix(
            "113576825867421328269022897092212996059028298184297174052451199438117344128029",
            10,
        )
        .unwrap();
        let z2: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let scalar1: BigInt = BigInt::from_str_radix(
            "51744903828951655645998112878951294312407757131136959404777653173288809286437",
            10,
        )
        .unwrap();
        let scalar2: BigInt = BigInt::from_str_radix(
            "47255149894527971878307606101838196302615590688404325145363490576551532104959",
            10,
        )
        .unwrap();

        let x_res: BigInt = BigInt::from_str_radix(
            "4553550298606723614503613575355345616769837055925158630918861136516274250754",
            10,
        )
        .unwrap();
        let y_res: BigInt = BigInt::from_str_radix(
            "55240462671849106909136195492796706817980992148608231216795717302256464992045",
            10,
        )
        .unwrap();
        let z_res: BigInt = BigInt::from_str_radix("1", 10).unwrap();

        let point1: ProjectivePoint = ProjectivePoint::new(x1.clone(), y1.clone(), z1.clone());
        let point2: ProjectivePoint = ProjectivePoint::new(x2.clone(), y2.clone(), z2.clone());

        let computed_result_projective_point: ProjectivePoint = curve.add(
            &curve.scalar_mul(&scalar1, &point1),
            &curve.scalar_mul(&scalar2, &point2),
        );

        let computed_result_affine_point: AffinePoint =
            computed_result_projective_point.to_affine(&curve).unwrap();

        let resut: ProjectivePoint = computed_result_affine_point.to_projective();

        assert_eq!(resut.x, x_res);
        assert_eq!(resut.y, y_res);
        assert_eq!(resut.z, z_res);
    }
}
