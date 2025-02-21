// helpers.rs
use super::*;
use crate::Curve;
use num_bigint::BigInt;
use num_traits::Zero;
use rand::{prelude::ThreadRng, Rng};

pub fn module(x: &BigInt, modulo: &BigInt) -> BigInt {
    let result: BigInt = x % modulo;

    if result < BigInt::zero() {
        return modulo + result;
    }

    result
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

pub fn get_curve() -> Curve {
    let p: BigInt = BigInt::from_str_radix(
        "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
        16,
    )
    .unwrap();
    let a: BigInt = BigInt::from_str_radix(
        "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
        16,
    )
    .unwrap();
    let b: BigInt = BigInt::from_str_radix(
        "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        16,
    )
    .unwrap();

    Curve::new(p, a, b).unwrap()
}

pub fn measure_average_execution_time<F>(method_name: &str, method_to_run: F, iterations: usize)
where
    F: Fn() -> u128,
{
    let mut total_duration_in_nano: u128 = 0u128;

    for _ in 0..iterations {
        total_duration_in_nano += method_to_run();
    }

    let avg_duration_in_nano: u128 = total_duration_in_nano / iterations as u128;
    let avg_duration_in_millis: u128 = avg_duration_in_nano / 1_000_000;

    println!();
    println!("Operation called: {}", method_name);
    println!("Average Execution time (ns): {}", avg_duration_in_nano);
    println!("Average Execution time (ms): {}", avg_duration_in_millis);
}
