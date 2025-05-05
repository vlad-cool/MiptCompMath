pub fn unflatten<const A: usize, const B: usize>(flat: &[f64]) -> [[f64; A]; B] {
    assert!(flat.len() == A * B);
    std::array::from_fn(|i| std::array::from_fn(|j| flat[i * A + j]))
}

pub fn close_enough(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

pub fn close_enough_arr<const N: usize>(a: &[f64; N], b: &[f64; N], epsilon: f64) -> bool {
    let mut sum: f64 = 0f64;
    for i in 0..N {
        sum += (a[i] - b[i]).abs();
    }

    close_enough(sum, 0f64, epsilon)
}

pub fn derivative<F>(f: &mut F, x: f64, h: f64) -> f64
where
    F: FnMut(f64) -> f64,
{
    (f(x + h) - f(x - h)) / (2f64 * h)
}

pub fn partial_derivative<F, const N: usize>(f: &mut F, x: &[f64; N], i: usize, j: usize, h: f64) -> f64
where
    F: FnMut(&[f64; N]) -> [f64; N],
{
    let mut x_1: [f64; N] = x.clone();
    let mut x_2: [f64; N] = x.clone();

    x_1[j] += h;
    x_2[j] -= h;

    (f(&x_1)[i] - f(&x_2)[i]) / (2f64 * h)
}

pub fn zero_enough_vec(a: &std::vec::Vec<f64>, epsilon: f64) -> bool {
    let mut sum: f64 = 0f64;
    for i in 0..a.len() {
        sum += a[i].abs();
    }

    close_enough(sum, 0f64, epsilon)
}

pub fn partial_derivative_vec<F>(f: &mut F, x: &std::vec::Vec<f64>, i: usize, j: usize, h: f64) -> f64
where
    F: FnMut(&std::vec::Vec<f64>) -> std::vec::Vec<f64>,
{
    let mut x_1: Vec<f64> = x.clone();
    let mut x_2: Vec<f64> = x.clone();

    x_1[j] += h;
    x_2[j] -= h;

    (f(&x_1)[i] - f(&x_2)[i]) / (2f64 * h)
}

pub fn c_n_k(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }

    let k = std::cmp::min(k, n - k);
    let mut res = 1;

    for i in 1..=k {
        res = res * (n - k + i) / i;
    }

    res
}

#[cfg(test)]
mod tests {
    use crate::utils::*;

    #[test]
    fn test_derivative() {
        fn f(x: f64) -> f64 {
            x.sin()
        }

        assert!(close_enough(derivative(&mut f, 0f64, 1e-6), 1f64, 1e-6));
        assert!(close_enough(
            derivative(&mut f, std::f64::consts::FRAC_PI_2, 1e-6),
            0f64,
            1e-6
        ));
    }

    #[test]
    fn test_partial_derivative() {
        fn f(x: &[f64; 2]) -> [f64; 2] {
            [x[0] * x[0] + x[1] * x[1], 2f64 * x[0] * x[1]]
        }

        assert!(close_enough(
            partial_derivative(&mut f, &[0f64, 0f64], 0, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[0f64, 0f64], 0, 1, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[0f64, 0f64], 1, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[0f64, 0f64], 1, 1, 1e-6),
            0f64,
            1e-6
        ));

        assert!(close_enough(
            partial_derivative(&mut f, &[1f64, 1f64], 0, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[1f64, 1f64], 0, 1, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[1f64, 1f64], 1, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(&mut f, &[1f64, 1f64], 1, 1, 1e-6),
            2f64,
            1e-6
        ));
    }

    #[test]
    fn test_c_n_k() {
        assert!(c_n_k(10, 0) == 1);
        assert!(c_n_k(10, 1) == 10);
        assert!(c_n_k(10, 2) == 45);
        assert!(c_n_k(10, 9) == 10);
        assert!(c_n_k(10, 10) == 1);
    }
}
