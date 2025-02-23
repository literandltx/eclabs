pub struct PrimeGenerator {
    current: u64,
}

#[allow(unused)]
impl PrimeGenerator {
    pub fn new() -> Self {
        PrimeGenerator { current: 1 }
    }

    pub fn next_prime(&mut self) -> u64 {
        self.current += 1;

        loop {
            if is_prime(&self.current) {
                let prime: u64 = self.current;
                self.current += 1;

                return prime;
            }

            self.current += 1;
        }
    }
}

// todo
#[allow(unused)]
fn is_prime(n: &u64) -> bool {
    if *n <= 1 {
        return false;
    }

    for i in 2..=(*n as f64).sqrt() as u64 {
        if n % i == 0 {
            return false;
        }
    }

    true
}
