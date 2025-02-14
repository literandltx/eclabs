// helpers.rs
use num_bigint::BigInt;
use num_traits::Zero;

pub fn module(x: &BigInt, modulo: &BigInt) -> BigInt {
    if *x < BigInt::zero() {
        return modulo + x;
    }

    x % modulo
}
