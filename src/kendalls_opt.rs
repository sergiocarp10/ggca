use std::cmp::Ordering;
use kendalls::Error;
use rayon::slice::ParallelSliceMut;

fn pairs_comparator(a: &f64, b: &f64) -> Ordering {
    a.partial_cmp(b).unwrap_or(Ordering::Greater)
}

pub fn tau_b_with_comparator(x: &[f64], y: &[f64]) -> Result<(f64, f64), Error> {

    if x.len() != y.len() {
        return Err(Error::DimensionMismatch {
            expected: x.len(),
            got: y.len(),
        });
    }

    if x.is_empty() {
        return Err(Error::InsufficientLength);
    }

    let n = x.len();

    let mut pairs: Vec<(f64, f64)> = x.iter().cloned().zip(y.iter().cloned()).collect();

    pairs.par_sort_by(|pair1, pair2| {
        let res = pairs_comparator(&pair1.0, &pair2.0);
        if res == Ordering::Equal {
            pairs_comparator(&pair1.1, &pair2.1)
        } else {
            res
        }
    });

    let mut v1_part_1 = 0usize;
    let mut v2_part_1 = 0isize;
    let mut tied_x_pairs = 0usize;
    let mut tied_xy_pairs = 0usize;
    let mut vt = 0usize;
    let mut consecutive_x_ties = 1usize;
    let mut consecutive_xy_ties = 1usize;

    for i in 1..n {
        let prev = &pairs[i - 1];
        let curr = &pairs[i];

        if curr.0 == prev.0 {
            consecutive_x_ties += 1;
            if curr.1 == prev.1 {
                consecutive_xy_ties += 1;
            } else {
                tied_xy_pairs += sum(consecutive_xy_ties - 1);
                consecutive_xy_ties = 1;
            }
        } else {
            update_x_group(
                &mut vt,
                &mut tied_x_pairs,
                &mut tied_xy_pairs,
                &mut v1_part_1,
                &mut v2_part_1,
                consecutive_x_ties,
                consecutive_xy_ties,
            );
            consecutive_x_ties = 1;
            consecutive_xy_ties = 1;
        }
    }

    update_x_group(
        &mut vt,
        &mut tied_x_pairs,
        &mut tied_xy_pairs,
        &mut v1_part_1,
        &mut v2_part_1,
        consecutive_x_ties,
        consecutive_xy_ties,
    );

    let mut swaps = 0usize;
    let mut pairs_dest: Vec<(f64, f64)> = vec![(0.0, 0.0); n];
    let mut segment_size = 1usize;

    while segment_size < n {
        for offset in (0..n).step_by(2 * segment_size) {
            let mut i = offset;
            let i_end = n.min(i + segment_size);
            let mut j = i_end;
            let j_end = n.min(j + segment_size);
            let mut copy_location = offset;

            while i < i_end || j < j_end {
                if i < i_end {
                    if j < j_end {
                        let a = &pairs[i].1;
                        let b = &pairs[j].1;
                        if pairs_comparator(a, b) == Ordering::Greater {
                            pairs_dest[copy_location] = pairs[j].clone();
                            j += 1;
                            swaps += i_end - i;
                        } else {
                            pairs_dest[copy_location] = pairs[i].clone();
                            i += 1;
                        }
                    } else {
                        pairs_dest[copy_location] = pairs[i].clone();
                        i += 1;
                    }
                } else {
                    pairs_dest[copy_location] = pairs[j].clone();
                    j += 1;
                }
                copy_location += 1;
            }
        }

        std::mem::swap(&mut pairs, &mut pairs_dest);

        segment_size <<= 1;
    }

    let mut v1_part_2 = 0usize;
    let mut v2_part_2 = 0isize;
    let mut tied_y_pairs = 0usize;
    let mut consecutive_y_ties = 1usize;
    let mut vu = 0usize;

    for j in 1..n {
        let prev = &pairs[j - 1];
        let curr = &pairs[j];

        if curr.1 == prev.1 {
            consecutive_y_ties += 1;
        } else {
            update_y_group(
                &mut vu,
                &mut tied_y_pairs,
                &mut v1_part_2,
                &mut v2_part_2,
                consecutive_y_ties,
            );
            consecutive_y_ties = 1;
        }
    }

    update_y_group(
        &mut vu,
        &mut tied_y_pairs,
        &mut v1_part_2,
        &mut v2_part_2,
        consecutive_y_ties,
    );

    // Generates T1 and T2 for significance
    let v1 = (v1_part_1 * v1_part_2) as f64;
    let v2 = (v2_part_1 * v2_part_2) as f64;

    // Prevents overflow on subtraction
    let num_pairs_f: f64 = ((n * (n - 1)) as f64) / 2.0; // sum(n - 1).as_();
    let tied_x_pairs_f: f64 = tied_x_pairs as f64;
    let tied_y_pairs_f: f64 = tied_y_pairs as f64;
    let tied_xy_pairs_f: f64 = tied_xy_pairs as f64;
    let swaps_f: f64 = (2 * swaps) as f64;

    let concordant_minus_discordant = num_pairs_f - tied_x_pairs_f - tied_y_pairs_f + tied_xy_pairs_f - swaps_f;
    let non_tied_pairs_multiplied = (num_pairs_f - tied_x_pairs_f) * (num_pairs_f - tied_y_pairs_f);
    let tau_b = concordant_minus_discordant / non_tied_pairs_multiplied.sqrt();

    // Significance
    let v0 = (n * (n - 1)) * (2 * n + 5);
    let n_f = n as f64;

    let v0_isize = v0 as isize;
    let vt_isize = vt as isize;
    let vu_isize = vu as isize;
    let var_s = (v0_isize - vt_isize - vu_isize) as f64 / 18.0
        + v1 / (2.0 * n_f * (n_f - 1.0))
        + v2 / (9.0 * n_f * (n_f - 1.0) * (n_f - 2.0));

    let s = tau_b * non_tied_pairs_multiplied.sqrt();
    let z = s / var_s.sqrt();

    // Limit range to fix computational errors
    Ok((tau_b.max(-1.0).min(1.0), z))

}

// Copied from Kendalls Lib

#[inline]
fn sum(n: usize) -> usize {
    n * (n + 1_usize) / 2_usize
}

/// Updated vt, v1_part_1, v2_part_1, tied_x_pairs, tied_xy_pairs variables with current tied group in X
fn update_x_group(
    vt: &mut usize,
    tied_x_pairs: &mut usize,
    tied_xy_pairs: &mut usize,
    v1_part_1: &mut usize,
    v2_part_1: &mut isize,
    consecutive_x_ties: usize,
    consecutive_xy_ties: usize,
) {
    *vt += consecutive_x_ties * (consecutive_x_ties - 1) * (2 * consecutive_x_ties + 5);
    *v1_part_1 += consecutive_x_ties * (consecutive_x_ties - 1);

    let consecutive_x_ties_i = consecutive_x_ties as isize;
    *v2_part_1 += consecutive_x_ties_i * (consecutive_x_ties_i - 1) * (consecutive_x_ties_i - 2);

    *tied_x_pairs += sum(consecutive_x_ties - 1);
    *tied_xy_pairs += sum(consecutive_xy_ties - 1);
}

/// Updated vu, tied_y_pairs, v1_part_2 and v2_part_2 variables with current tied group in Y
fn update_y_group(
    vu: &mut usize,
    tied_y_pairs: &mut usize,
    v1_part_2: &mut usize,
    v2_part_2: &mut isize,
    consecutive_y_ties: usize,
) {
    *vu += consecutive_y_ties * (consecutive_y_ties - 1) * (2 * consecutive_y_ties + 5);
    *v1_part_2 += consecutive_y_ties * (consecutive_y_ties - 1);

    let consecutive_y_ties_i = consecutive_y_ties as isize;
    *v2_part_2 += consecutive_y_ties_i * (consecutive_y_ties_i - 1) * (consecutive_y_ties_i - 2);

    *tied_y_pairs += sum(consecutive_y_ties - 1);
}