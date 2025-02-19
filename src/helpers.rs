// helpers.rs
use num_bigint::BigInt;
use num_traits::Zero;
use rand::{prelude::ThreadRng, Rng};
use std::time::Instant;

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

fn measure_average_execution_time<F>(method_name: &str, method_to_run: F, iterations: usize)
where
    F: Fn() -> (),
{
    let mut total_duration_in_nano: u128 = 0u128;

    for _ in 0..iterations {
        let start_time: Instant = Instant::now();

        method_to_run();

        let end_time: Instant = Instant::now();
        total_duration_in_nano += end_time.duration_since(start_time).as_nanos();
    }

    let avg_duration_in_nano: u128 = total_duration_in_nano / iterations as u128;
    let avg_duration_in_millis: u128 = avg_duration_in_nano / 1_000_000;

    println!("Operation called: {}", method_name);
    println!("Average Execution time (ns): {}", avg_duration_in_nano);
    println!("Average Execution time (ms): {}", avg_duration_in_millis);
    println!();
}
