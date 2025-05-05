pub fn tensor_product(
    l: &std::vec::Vec<std::vec::Vec<f64>>,
    r: &std::vec::Vec<std::vec::Vec<f64>>,
) -> std::vec::Vec<std::vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let (l, r) = (r, l);

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len() * r[0].len()]; l.len() * r.len()];

    for a in 0..r.len() {
        for b in 0..r[0].len() {
            for c in 0..l.len() {
                for d in 0..l[0].len() {
                    res[a * l.len() + c][b * l[0].len() + d] = l[c][d] * r[a][b];
                }
            }
        }
    }

    res
}

pub fn matrix_product(
    l: &std::vec::Vec<std::vec::Vec<f64>>,
    r: &std::vec::Vec<std::vec::Vec<f64>>,
) -> std::vec::Vec<std::vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    assert!(l[0].len() == r.len());

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; r[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..r[0].len() {
            for c in 0..r.len() {
                res[a][b] += r[c][b] * l[a][c];
            }
        }
    }

    res
}

pub fn matrix_sum(l: &std::vec::Vec<std::vec::Vec<f64>>, r: &std::vec::Vec<std::vec::Vec<f64>>) -> std::vec::Vec<std::vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    assert!(l.len() == r.len());
    assert!(l[0].len() == r[0].len());
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..l[0].len() {
            res[a][b] = l[a][b] + r[a][b];
        }
    }

    res
}

pub fn matrix_sub(l: &std::vec::Vec<std::vec::Vec<f64>>, r: &std::vec::Vec<std::vec::Vec<f64>>) -> std::vec::Vec<std::vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    assert!(l.len() == r.len());
    assert!(l[0].len() == r[0].len());
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..l[0].len() {
            res[a][b] = l[a][b] - r[a][b];
        }
    }

    res
}

pub fn matrix_transpose(a: &std::vec::Vec<std::vec::Vec<f64>>) -> std::vec::Vec<std::vec::Vec<f64>> {
    for i in 0..(a.len() - 1) {
        assert!(a[i].len() == a[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; a.len()]; a[0].len()];

    for i in 0..a.len() {
        for j in 0..a[0].len() {
            res[j][i] = a[i][j];
        }
    }

    res
}

pub fn matrix_scale(l: &std::vec::Vec<std::vec::Vec<f64>>, r: f64) -> std::vec::Vec<std::vec::Vec<f64>> {
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for i in 0..l.len() {
        for j in 0..l[0].len() {
            res[i][j] = l[i][j] * r;
        }
    }

    res
}

/// Tests

#[cfg(test)]
mod tests {
    use crate::matrix_operations::*;

    #[test]
    fn test_tensor_product() {
        let a: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0],
            vec![2.0, 0.0, 4.0],
            vec![3.0, 3.0, 0.0],
        ];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

        let result: Vec<Vec<f64>> = tensor_product(&b, &a);

        let reference: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
            vec![2.0, 0.0, 4.0, 4.0, 0.0, 8.0],
            vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
            vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
            vec![3.0, 0.0, 6.0, 4.0, 0.0, 8.0],
            vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
        ];

        println!("{:?}", result);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 3.0, 0.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

        let result: Vec<Vec<f64>> = tensor_product(&b, &a);

        let reference: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
            vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
            vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
            vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
        ];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_product() {
        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 2.0, 0.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 1.0], vec![1.5, 2.0], vec![2.0, 1.0]];

        let result: Vec<Vec<f64>> = matrix_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![vec![1.5, 2.0], vec![6.0, 7.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![1.0, 1.0, 1.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0], vec![1.5], vec![2.0]];

        let result: Vec<Vec<f64>> = matrix_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![vec![1.5], vec![4.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_sum_sub() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let b: Vec<Vec<f64>> = vec![vec![0.5, 1.5, 2.5], vec![3.5, 4.5, 3.5]];

        let result: Vec<Vec<f64>> = matrix_sum(&a, &b);
        let reference: Vec<Vec<f64>> = vec![vec![1.5, 3.5, 5.5], vec![7.5, 9.5, 9.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let result: Vec<Vec<f64>> = matrix_sub(&a, &b);
        let reference: Vec<Vec<f64>> = vec![vec![0.5, 0.5, 0.5], vec![0.5, 0.5, 2.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_transpose() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let result: Vec<Vec<f64>> = matrix_transpose(&a);
        let reference: Vec<Vec<f64>> = vec![vec![1.0, 4.0], vec![2.0, 5.0], vec![3.0, 6.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_scale_transpose() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let result: Vec<Vec<f64>> = matrix_scale(&a, 0.5);
        let reference: Vec<Vec<f64>> = vec![vec![0.5, 1.0, 1.5], vec![2.0, 2.5, 3.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

}
