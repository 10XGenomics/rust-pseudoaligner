use std::collections::HashMap;
use bam::EqClass;

pub struct EqClassCounts {
    pub(crate) counts: HashMap<EqClass, u32>,
}

impl EqClassCounts {
    pub fn new() -> EqClassCounts {
        EqClassCounts {
            counts: HashMap::new(),
        }
    }
}

pub fn em(nitems: usize, eqclasses: &EqClassCounts) -> Vec<f64> {

    // initialize weights
    let mut weights = vec![1.0/(nitems as f64); nitems];
    
    let mut iters = 0;

    // Abundance required to 'care' about a relative change
    let rel_care_thresh = 1e-3 / (nitems as f64);

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

        for i in 0 .. nitems {
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
            simpsons += new_weights*new_weights;
        }

        println!("iter: {}, div: {}, rel_diff: {}, abs_diff: {}", iters, 1.0/simpsons, max_rel_diff, max_abs_diff);
        iters += 1;
        if (max_abs_diff < 0.0001 && max_rel_diff < 0.00025) || iters > 5000 {
            break;
        }

        
    }

    weights
}




#[cfg(test)]
mod test {
    use super::*;
    use bam::EqClass;

    #[test]
    fn simple_inf() {

        let mut counts = HashMap::new();

        let eqA = EqClass::from(vec![0]);
        let eqAB = EqClass::from(vec![0,1]);

        let eqC = EqClass::from(vec![2]);
        let eqD = EqClass::from(vec![3]);

        counts.insert(eqA, 1);
        counts.insert(eqAB, 19);
        counts.insert(eqC, 10);
        counts.insert(eqD, 10);

        let eqc = EqClassCounts { counts };

        let res = em(4, &eqc);

        println!("{:?}", res);
    }
}