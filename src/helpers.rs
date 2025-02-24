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
    let num_bytes: usize = (bit_length + 7) / 8;
    let mask: BigInt = (BigInt::one() << bit_length) - BigInt::one();

    let mut rng: ThreadRng = rand::rng();
    let mut result: BigInt = BigInt::zero();

    let random_bytes: Vec<u8> = (0..num_bytes).map(|_| rng.random()).collect();

    for (i, byte) in random_bytes.iter().enumerate() {
        result += BigInt::from(*byte) << (i * 8);
    }

    result & mask
}

pub fn get_curve_p256() -> Curve {
    let p: BigInt = BigInt::from_str_radix(
        "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
        16,
    )
    .expect("Failed to parse curve parameter p");
    let a: BigInt = BigInt::from_str_radix(
        "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
        16,
    )
    .expect("Failed to parse curve parameter a");
    let b: BigInt = BigInt::from_str_radix(
        "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
        16,
    )
    .expect("Failed to parse curve parameter b");

    Curve::new(p, a, b)
        .expect("Failed to create curve")
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
    let avg_duration_in_millis: f64 = avg_duration_in_nano as f64 / 1_000_000.;

    println!();
    println!("Operation called: {}", method_name);
    println!("Average Execution time (ns): {}", avg_duration_in_nano);
    println!("Average Execution time (ms): {}", avg_duration_in_millis);
}
