// helpers.rs
use num_bigint::BigInt;
use num_traits::Zero;
use rand::prelude::ThreadRng;
use rand::Rng;

pub fn module(x: &BigInt, modulo: &BigInt) -> BigInt {
    if *x < BigInt::zero() {
        return modulo + x;
    }

    x % modulo
}

pub fn generate_random_bigint(bit_length: usize) -> BigInt {
    let mut rng: ThreadRng = rand::rng();
    let mut result: BigInt = BigInt::zero();

    for i in 0..bit_length {
        let bit: bool = rng.random();

        if bit {
            result += BigInt::from(1u32) << i;
        }
    }

    result
}
