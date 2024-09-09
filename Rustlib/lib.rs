use pyo3::prelude::*;
use std::collections::HashMap;
#[pyfunction]
fn calculate_distance_matrix(molecules: Vec<Vec<f64>>, threshold: f64)->Vec<(usize, usize, f64)>{
    let n = molecules.len();
    let t = molecules[0].len();
    let mut distance_full = vec![vec![0.0;n];n];
    let mut distance_sums: HashMap<(usize, usize), f64> = HashMap::new();
// Previously I calculated all distances but that was slow af. Now we only calculate the summed distance
// between each pair of molecules. This is intended to be quicker.
    for i in 0..n{
        for j in i+1..n{
            let mut dist_sum = 0.0;
            for column in 0..t{
                let dist = (molecules[i][column] - molecules[j][column]).powi(2).sqrt();
                dist_sum +=dist;
                distance_full[i][j] = dist;
                distance_full[j][i] = dist;
            }
            distance_sums.insert((i,j),dist_sum);
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
        .collect::<Vec<(usize, usize, f64)>>()
}

#[pyfunction]
fn dynamic_time_warping(molecules: Vec<Vec<f64>>, threshold: f64) -> Vec<(Vec<(usize, usize)>, f64)>{
    let n = molecules.len();
    let t = molecules[0].len();
    let mut cost_map = Vec::new();

    for i in 0..n{
        for j in i+1..n{
            let mut costmatrix = vec![vec![f64::INFINITY; t+1]; t+1];
            costmatrix[0][0] = 0.0;

            for k in 0..t{
                for l in 0..t{
                    let cost = molecules[i][k].abs() - molecules[j][l].abs();
                    costmatrix[k+1][l+1]= cost + costmatrix[k][l].min(costmatrix[k+1][l].min(costmatrix[k][l+1]));
                }
            }
            // Backtrack to determine optimal path
            let mut matches = Vec::new();
            let (mut k, mut l) = (t-1, t-1);

            while k>0 || l>0 {
                matches.push((k,l));
                let option_diag = if k > 0 && l > 0 {costmatrix[k][l]} else {f64::INFINITY};
                let option_up = if k > 0 {costmatrix[k][l+1]} else {f64::INFINITY};
                let option_left = if l > 0 {costmatrix[k+1][l]} else {f64::INFINITY};
                let best_move = vec![option_diag, option_up, option_left].into_iter().enumerate().min_by(|a,b| a.1.partial_cmp(&b.1).unwrap()).unwrap().0;
                match best_move{
                    0=> {k-=1; l -=1;}
                    1=> {k -= 1;}
                    2=> {l -=1;}
                    _=> {}
                }
            }
            matches.reverse();

            let final_cost = costmatrix[t][t];
            if final_cost <= threshold{
                cost_map.push((matches, final_cost));
            }
        }
    }
    cost_map
}
#[pymodule]
fn distance_calc(_py: Python, m: &PyModule)-> PyResult<()>{
    m.add_function(wrap_pyfunction!(calculate_distance_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(dynamic_time_warping, m)?)?;
    Ok(())
}