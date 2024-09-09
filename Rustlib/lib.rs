use pyo3::prelude::*;
use std::collections::HashMap;
#[pyfunction]
fn calculate_distance_matrix(molecules: Vec<Vec<f64>>, threshold: f64)->Vec<(usize, usize, f64)>{
    let n = molecules.len();
    let t = molecules[0].len();
    let mut distance_sums: HashMap<(usize, usize), f64> = Hashmap::new();
    let mut idx = 0;
// Previously I calculated all distances but that was slow af. Now we only calculate the summed distance
// between each pair of molecules. This is intended to be quicker.
    for column in 0..t{
        for i in 0..n{
            for j in i + 1..n{
                let dist = (molecules[i][column] - molecules[j][column]).powi(2).sqrt();
                let pair = (i, j);
                *distance_sums.entry(pair).or_insert(0.0) += dist;
            }
        }
    }
    let mut average_distances: Vec<(usize, usize, f64)> = Vec::new();
    for ((i,j), dist_sum) in distance_sums.iter(){
    let avg_dist = dist_sum/ t as f64;
    average_distances.push((*i,*j, avg_dist));
    }
    // Made another variable so av_dist is not sorted. Is this necessary? Would probs be faster without extra variable.
    let mut sorted_distances = average_distances.clone();
    sorted_distances.sort_by(|a,b| a.2.partial_cmp(&b.2).unwrap());
    let threshold_index = (threshold * sorted_distances.len() as f64).floor() as usize;
    let threshold_value = sorted_distances[threshold_index].2;
    average_distances
        .into_iter()
        .filter(|&(_,_,avg_dist)|avg_dist <= threshold_value)
        .collect()
}
#[pymodule]
fn distance_calc(_py: Python, m: &PyModule)-> PyResult<()>{
    m.add_function(wrap_pyfunction!(calculate_distance_matrix, m)?)?;
    Ok(())
}