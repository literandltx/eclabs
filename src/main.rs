mod helpers;
mod p256_test;
mod prime_generator;

use crate::helpers::{generate_random_bigint, module};
use num::Integer;
use num_bigint::BigInt;
use num_bigint::RandBigInt;
use num_traits::{Num, One, Pow, Zero};
use rand::Rng;
use rs_sha256;
use std::hash::BuildHasher;
use std::hash::{Hash, Hasher};

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

#[derive(Clone, Debug)]
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
    fn new(p: BigInt, a: BigInt, b: BigInt) -> Option<Curve> {
        let d: BigInt = a.clone().pow(3u32) * 4 + b.clone().pow(2u32) * 27u32;

        match d != BigInt::zero() {
            true => Some(Curve { p, a, b }),
            false => None,
        }
    }

    fn verify_affine_point(&self, point: &AffinePoint) -> bool {
        let left: BigInt = point.clone().y.pow(2u32);
        let right: BigInt = point.clone().x.pow(3u32) + &self.a * point.clone().x + &self.b;

        left % &self.p == right % &self.p
    }

    fn verify_projective_point(&self, point: &ProjectivePoint) -> bool {
        let z1: BigInt = point.z.clone();
        let z2: BigInt = point.z.clone().pow(2u32);
        let z3: BigInt = point.z.clone().pow(3u32);

        let left: BigInt = point.clone().y.pow(2u32) * z1;
        let right: BigInt =
            point.clone().x.pow(3u32) + &self.a * z2 * point.clone().x + &self.b * z3;

        left % &self.p == right % &self.p
    }

    fn add(&self, a: &ProjectivePoint, b: &ProjectivePoint) -> ProjectivePoint {
        if a.x == BigInt::zero() && a.y == BigInt::one() && a.z == BigInt::zero() {
            return b.clone();
        }

        if b.x == BigInt::zero() && b.y == BigInt::one() && b.z == BigInt::zero() {
            return a.clone();
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
        let x: BigInt = &v * &_a;

        // Y3 := U*(V^2*V2 - A) - V^3*U2
        let y: BigInt = &u * (&v * &v * &v2 - &_a) - &v * &v * &v * &u2;

        // Z3 := V^3*W
        let z: BigInt = &v * &v * &v * &w;

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
        let x: BigInt = 2 * &h * &s;

        // Y' := W*(4*B - H) - 8*Y^2*S^2
        let y: BigInt = &w * (4 * &b - &h) - 8 * &point.y * &point.y * &s * &s;

        // Z' := 8*S^3
        let z: BigInt = 8 * &s * &s * &s;

        ProjectivePoint::new(
            module(&x, &self.p),
            module(&y, &self.p),
            module(&z, &self.p),
        )
    }

    fn scalar_mul(&self, scalar: &BigInt, point: &ProjectivePoint) -> ProjectivePoint {
        let mut result: ProjectivePoint = ProjectivePoint::neutral();
        let mut temp: ProjectivePoint = point.clone();

        for i in scalar.to_str_radix(2).into_bytes().iter().rev() {
            if *i - b'0' == 1 {
                result = self.add(&result, &temp);
            }

            temp = self.double(&temp);
        }

        result
    }

    fn scalar_mul_montgomery(&self, scalar: &BigInt, point: &ProjectivePoint) -> ProjectivePoint {
        let mut r0: ProjectivePoint = ProjectivePoint::neutral();
        let mut r1: ProjectivePoint = point.clone();

        for i in scalar.to_str_radix(2).into_bytes().iter() {
            if *i - b'0' == 0 {
                r1 = self.add(&r1, &r0);
                r0 = self.double(&r0);
            } else {
                r0 = self.add(&r0, &r1);
                r1 = self.double(&r1);
            }
        }

        r0
    }

    fn get_generator(&self) -> ProjectivePoint {
        let gen_x: BigInt = BigInt::from_str_radix(
            "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
            16,
        )
        .unwrap();
        let gen_y: BigInt = BigInt::from_str_radix(
            "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5",
            16,
        )
        .unwrap();

        ProjectivePoint::new(gen_x, gen_y, BigInt::one())
    }

    fn get_order(&self) -> BigInt {
        BigInt::from_str_radix(
            "115792089210356248762697446949407573529996955224135760342422259061068512044369",
            10,
        )
        .unwrap()
    }

    fn generate_random_projective_point_p256_naive(&self) -> ProjectivePoint {
        let rand: BigInt = generate_random_bigint(128);

        assert!(
            self.verify_projective_point(&self.get_generator()),
            "The generator do not lay on current curve"
        );

        let random_point: ProjectivePoint = self
            .scalar_mul_montgomery(&rand, &self.get_generator())
            .to_affine(&self)
            .unwrap()
            .to_projective();

        assert!(
            self.verify_projective_point(&random_point),
            "The point do not lay on current curve"
        );

        random_point
    }

    fn generate_random_projective_point_p256_fast(&self) -> ProjectivePoint {
        loop {
            let x: BigInt = generate_random_bigint(256);
            let right: BigInt = x.clone().pow(3u32) + &self.a * x.clone() + &self.b;

            // right^((p - 1) / 2)
            let check_value: BigInt =
                right.modpow(&((&self.p - BigInt::one()) / BigInt::from(2u32)), &self.p);

            // check is there existing solutions
            if check_value.eq(&BigInt::one()) {
                let k: BigInt = (&self.p - BigInt::from(3u32)) / BigInt::from(4u32);

                // y = + right^(k + 1)
                let y: BigInt = right.modpow(&(k + BigInt::one()), &self.p);

                return ProjectivePoint::new(x, y, BigInt::one());
            }
        }
    }
}

struct Actor {
    curve: Curve,
    sk: BigInt,
}

impl Actor {
    fn new(curve: Curve) -> Actor {
        let mut rng = rand::thread_rng();
        let sk: BigInt =
            RandBigInt::gen_bigint_range(&mut rng, &BigInt::from(2u32), &curve.get_order());
        Actor { curve, sk }
    }

    fn generate_secret_key(&self) -> BigInt {
        let mut rng = rand::thread_rng();
        RandBigInt::gen_bigint_range(&mut rng, &BigInt::from(2u32), &self.curve.get_order())
    }

    fn get_public_key(&self, base_point: &ProjectivePoint) -> ProjectivePoint {
        self.curve.scalar_mul(&self.sk, base_point)
    }

    fn update_secret_key(&mut self, base_point: &ProjectivePoint) -> ProjectivePoint {
        self.sk = self.generate_secret_key();
        self.curve.scalar_mul(&self.sk, &base_point)
    }

    fn compute_common_secret(&self, key: &BigInt, pre_key: &ProjectivePoint) -> BigInt {
        let common_point = self.curve.scalar_mul(key, &pre_key).to_affine(&self.curve);

        common_point.unwrap().x
    }

    fn encrypt_aes(key_str: String, plaintext: String) -> Vec<u8> {
        let key_bytes = key_str.into_bytes();
        let plaintext_bytes = plaintext.into_bytes();

        let res: Vec<u8> = key_bytes
            .iter()
            .zip(plaintext_bytes.iter())
            .map(|(&x1, &x2)| x1 ^ x2)
            .collect();

        res
    }

    fn decrypt_aes(key_str: String, encrypted_data: Vec<u8>) -> String {
        let key_bytes = key_str.into_bytes();

        let plaintext: Vec<u8> = key_bytes
            .iter()
            .zip(encrypted_data.iter())
            .map(|(&x1, &x2)| x1 ^ x2)
            .collect();

        String::from_utf8(plaintext).expect("failed to convert vector of bytes to string")
    }

    fn encapsulate(key_str: String, plaintext: String) -> Vec<u8> {
        Actor::encrypt_aes(key_str, plaintext)
    }

    fn decapsulate(key_str: String, encrypted_data: Vec<u8>) -> String {
        Actor::decrypt_aes(key_str, encrypted_data)
    }

    fn hash(message: String) -> u64 {
        let mut sha256hasher = rs_sha256::Sha256State::default().build_hasher();
        sha256hasher.write(message.as_bytes());
        sha256hasher.finish()
    }

    fn encrypt(
        &self,
        message: String,
        receiver_public_key: &ProjectivePoint,
        base_point: &ProjectivePoint,
    ) -> (ProjectivePoint, Vec<u8>, Vec<u8>) {
        let k = self.generate_secret_key();
        let encrypted_m = Actor::encrypt_aes(k.to_string(), message);

        let ephemeral_key: BigInt = self.generate_secret_key();
        let public_ephemeral_key = self.curve.scalar_mul(&ephemeral_key, &base_point);
        let common_key = self.compute_common_secret(&ephemeral_key, &receiver_public_key);
        let encapsulated_key = Actor::encapsulate(common_key.to_string(), k.to_string());
        (public_ephemeral_key, encapsulated_key, encrypted_m)
    }

    fn decrypt(&self, ciphertext: (ProjectivePoint, Vec<u8>, Vec<u8>)) -> String {
        let (public_ephemeral_key, encapsulated_key, encrypted_m): (
            ProjectivePoint,
            Vec<u8>,
            Vec<u8>,
        ) = ciphertext;

        let common_key: BigInt = self.compute_common_secret(&self.sk, &public_ephemeral_key);

        let k = Actor::decapsulate(common_key.to_string(), encapsulated_key);

        Actor::decrypt_aes(k.to_string(), encrypted_m)
    }

    fn sign(&self, message: String, base_point: &ProjectivePoint) -> (BigInt, BigInt) {
        let h = Actor::hash(message.to_string());
        let n = self.curve.get_order();

        // k, x are pre-initialised to avoid error, they will be changed anyway
        let mut k: BigInt = BigInt::zero();
        let mut x: BigInt = BigInt::zero();
        let mut r: BigInt = BigInt::zero();
        while r.eq(&BigInt::ZERO) {
            k = self.generate_secret_key();
            x = self
                .curve
                .scalar_mul(&k, base_point)
                .to_affine(&self.curve)
                .unwrap()
                .x;
            (_, r) = x.div_mod_floor(&n);
        }

        let (_, s) =
            (BigInt::modinv(&k, &n).unwrap() * (h + self.sk.clone() * r.clone())).div_mod_floor(&n);

        (r, s)
    }

    fn verify(
        &self,
        sign: &(BigInt, BigInt),
        message: String,
        signer_pk: &ProjectivePoint,
        base_point: &ProjectivePoint,
    ) -> bool {
        let h: u64 = Actor::hash(message.to_string());
        let n: BigInt = self.curve.get_order();
        let (r, s): &(BigInt, BigInt) = sign;

        let (_, u1): (BigInt, BigInt) =
            (BigInt::modinv(&s, &n).unwrap() * h.clone()).div_mod_floor(&n);
        let (_, u2): (BigInt, BigInt) =
            (BigInt::modinv(&s, &n).unwrap() * r.clone()).div_mod_floor(&n);

        let (_, v) = self
            .curve
            .add(
                &self.curve.scalar_mul(&u1, base_point),
                &self.curve.scalar_mul(&u2, signer_pk),
            )
            .to_affine(&self.curve)
            .unwrap()
            .x
            .div_mod_floor(&n);

        v.eq(r)
    }
}

fn main() {}
