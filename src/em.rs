use crate::bam::EqClass;
use std::collections::HashMap;

pub struct EqClassCounts {
    pub nitems: usize,
    pub counts: HashMap<EqClass, u32>,
}

impl EqClassCounts {
    pub fn new() -> EqClassCounts {
        EqClassCounts {
            nitems: 0,
            counts: HashMap::new(),
        }
    }
}

impl EmProblem for EqClassCounts {
    fn init(&self) -> Vec<f64> {
        vec![1.0 / (self.nitems as f64); self.nitems]
    }

    fn reg(&self, theta: &mut [f64]) {
        let mut norm = 0.0;

        for i in 0..theta.len() {
            // Clamp weight
            let mut v = theta[i];
            if v > 1.0 {
                v = 1.0
            };
            if v < 0.0 {
                v = 1e-15
            };
            theta[i] = v;
            norm += v;
        }

        let inv_norm = 1.0 / norm;
        for i in 0..theta.len() {
            theta[i] = theta[i] * inv_norm;
        }
    }

    fn F(&self, theta1: &[f64], theta2: &mut [f64]) {
        let nitems = self.nitems;

        let mut total_counts = 0.0;

        for i in 0..theta2.len() {
            theta2[i] = 0.0;
        }

        for (class, count) in &self.counts {
            let mut norm = 0.0;
            for tx in class {
                norm += theta1[*tx as usize];
            }

            for tx in class {
                let tx_count = theta1[*tx as usize] / norm * (*count as f64);
                theta2[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;

        for i in 0..nitems {
            let old_weights = theta1[i];
            let new_weights = theta2[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            theta2[i] = new_weights;
        }
    }

    fn L(&self, theta: &[f64]) -> f64 {
        let mut ll = 0.0;

        for (class, count) in &self.counts {
            // Exclude empty equivalence classes
            if class.len() == 0 {
                continue;
            }

            let mut theta_sum: f64 = 0.0;
            for tx in class {
                theta_sum += theta[*tx as usize];
            }

            ll += (*count as f64) * theta_sum.ln();
        }

        ll
    }
}

pub fn em(eqclasses: &EqClassCounts) -> Vec<f64> {
    let nitems = eqclasses.nitems;

    // initialize weights
    let mut weights = vec![1.0 / (nitems as f64); nitems];

    let mut iters = 0;

    // Abundance required to 'care' about a relative change
    let _rel_care_thresh = 1e-3 / (nitems as f64);

    loop {
        let mut pseudocounts = vec![0.0; nitems];
        let mut total_counts = 0.0;

        for (class, count) in &eqclasses.counts {
            let mut norm = 0.0;
            for tx in class {
                norm += weights[*tx as usize];
            }

            for tx in class {
                let tx_count = weights[*tx as usize] / norm * (*count as f64);
                pseudocounts[*tx as usize] += tx_count;
                total_counts += tx_count;
            }
        }

        let mut max_abs_diff = 0.0;
        let mut max_rel_diff = 0.0;
        let mut simpsons = 0.0;

        for i in 0..nitems {
            let old_weights = weights[i];
            let new_weights = pseudocounts[i] / total_counts;

            let abs_diff = (old_weights - new_weights).abs();
            let rel_diff = abs_diff / old_weights;

            if abs_diff > max_abs_diff {
                max_abs_diff = abs_diff;
            }

            if new_weights > 1e-2 && rel_diff > max_rel_diff {
                max_rel_diff = rel_diff
            }

            weights[i] = new_weights;
            simpsons += new_weights * new_weights;
        }

        let ll = eqclasses.L(&weights);
        println!(
            "iter: {}, ll: {}, div: {}, rel_diff: {}, abs_diff: {}",
            iters,
            ll,
            1.0 / simpsons,
            max_rel_diff,
            max_abs_diff
        );
        iters += 1;
        if (max_abs_diff < 0.00005 && max_rel_diff < 0.0001) || iters > 5000 {
            break;
        }
    }

    weights
}

/// Encapsulate an EM optimization problem so that it can run through an accelerated EM loop (SquareM).
pub trait EmProblem {
    // Make an initial estimate of the parameter vector. May be a naive estimate.
    fn init(&self) -> Vec<f64>;

    // Regularize a parameter vector -- fixup an inconsistencies in the parameters.
    // E.g. legalize values or enforce global constrains.
    fn reg(&self, theta: &mut [f64]);

    // Update the parameters -- one EM step
    fn F(&self, theta1: &[f64], theta_2: &mut [f64]);

    // Compute the likelihood of a parameter set
    fn L(&self, theta: &[f64]) -> f64;
}

/// SquareM EM acceleration method.
/// As described in:
/// Varadhan, Ravi, and Christophe Roland.
/// "Simple and globally convergent methods for accelerating the convergence of any EM algorithm."
/// Scandinavian Journal of Statistics 35.2 (2008): 335-353.
/// Takes an implementation of `EmProblem` and applies the accelerated EM algorithm.
pub fn squarem<T: EmProblem>(p: &T) -> Vec<f64> {
    // Array for holding theta
    let mut theta0 = p.init();
    let mut theta1 = theta0.clone();
    let mut theta2 = theta0.clone();
    let mut theta_sq = theta0.clone();

    let mut r = theta0.clone();
    let mut v = theta0.clone();

    let n = theta0.len();
    let _prevLL = std::f64::MIN;
    let mut iters = 0;

    let rel_care_thresh = Some(1e-5);

    loop {
        // Get theta1
        p.F(&theta0, &mut theta1);
        p.F(&theta1, &mut theta2);

        let mut rsq: f64 = 0.0;
        let mut vsq: f64 = 0.0;

        for i in 0..n {
            r[i] = theta1[i] - theta0[i];
            rsq += r[i].powi(2);

            v[i] = theta2[i] - theta1[i] - r[i];
            vsq += v[i].powi(2);
        }

        let mut alpha = -rsq.sqrt() / vsq.sqrt();
        let mut alpha_tries = 1;
        let mut lsq = 0.0;
        let mut l2 = 0.0;

        loop {
            let alpha_sq = alpha.powi(2);

            for i in 0..n {
                theta_sq[i] = theta0[i] - 2.0 * alpha * r[i] + alpha_sq * v[i]
            }

            p.reg(&mut theta_sq);

            lsq = p.L(&theta_sq);
            l2 = p.L(&theta2);

            if lsq > l2 || alpha_tries > 5 {
                break;
            } else {
                alpha = (alpha + -1.0) / 2.0;
            }

            alpha_tries += 1;
        }

        let (max_rel_diff, max_abs_diff) = if lsq > l2 {
            let diff = diffs(&theta0, &theta_sq, rel_care_thresh);
            std::mem::swap(&mut theta0, &mut theta_sq);
            diff
        } else {
            let diff = diffs(&theta0, &theta2, rel_care_thresh);
            std::mem::swap(&mut theta0, &mut theta2);
            diff
        };

        println!(
            "iter: {}, ll2: {}, llsq: {}, alpha_tries: {}, rel_diff: {}, abs_diff: {}",
            iters, l2, lsq, alpha_tries, max_rel_diff, max_abs_diff
        );
        iters += 1;

        if (max_abs_diff < 5e-4 && max_rel_diff < 5e-3) || iters > 2000 {
            break;
        }
    }

    theta0
}

/// Compute the change in the parameter vectors, returning the largest relative and absolute change, respectively.
/// Only parameters with a value greater than rel_thresh (if set), are counted in the relative change check.
fn diffs(t1: &[f64], t2: &[f64], rel_thresh: Option<f64>) -> (f64, f64) {
    let mut max_abs_diff = 0.0;
    let mut max_rel_diff = 0.0;

    for i in 0..t1.len() {
        let old_weights = t1[i];
        let new_weights = t2[i];

        let abs_diff = (old_weights - new_weights).abs();
        let rel_diff = abs_diff / old_weights;

        if abs_diff > max_abs_diff {
            max_abs_diff = abs_diff;
        }

        if rel_thresh.map_or(true, |thresh| new_weights > thresh) && rel_diff > max_rel_diff {
            max_rel_diff = rel_diff
        }
    }

    (max_rel_diff, max_abs_diff)
}

#[cfg(test)]
mod test_em {
    use super::*;
    use crate::bam::EqClass;

    fn test_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eqA = EqClass::from(vec![0]);
        let eqAB = EqClass::from(vec![0, 1]);

        let eqC = EqClass::from(vec![2]);
        let eqD = EqClass::from(vec![3]);

        counts.insert(eqA, 1);
        counts.insert(eqAB, 19);
        counts.insert(eqC, 10);
        counts.insert(eqD, 10);

        EqClassCounts { counts, nitems: 4 }
    }

    fn test2_ds() -> EqClassCounts {
        let mut counts = HashMap::new();

        let eqA = EqClass::from(vec![0]);
        let eqAB = EqClass::from(vec![0, 1]);

        let eqC = EqClass::from(vec![2]);
        let eqD = EqClass::from(vec![3]);

        let eqE = EqClass::from(vec![4, 5]);

        counts.insert(eqA, 1);
        counts.insert(eqAB, 19);
        counts.insert(eqC, 10);
        counts.insert(eqD, 10);
        counts.insert(eqE, 20);

        EqClassCounts { counts, nitems: 6 }
    }

    #[test]
    fn simple_inf() {
        let eqc = test_ds();
        let res = em(&eqc);

        println!("{:?}", res);
    }

    #[test]
    fn med_inf() {
        let eqc = test2_ds();
        let res = em(&eqc);

        println!("{:?}", res);
    }

    #[test]
    fn accel_inf() {
        let eqc = test_ds();
        let res = squarem(&eqc);

        println!("{:?}", res);
    }
}
