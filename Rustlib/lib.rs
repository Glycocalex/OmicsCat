use pyo3::prelude::*;
#[pyfunction]
fn calculate_distance_matrix(molecules: Vec<Vec<f64>>)->Vec<f64>{
    let n = molecules.len();
    let mut distances = Vec::with_capacity(n*(n-1)/2);
    let mut idx = 0;

    for i in 0..n{
        for j in i + 1..n{
            let dist = molecules[i]
                .iter()
                .zip(molecules[j].iter())
                .map(|(a,b)| (a-b).powi(2))
                .sum::<f64>()
                .sqrt();
            distances[idx] = dist;
            idx +=1;
        }
    }
    distances
}
#[pymodule]
fn distance_calc(_py: Python, m: &PyModule)-> PyResult<()>{
    m.add_function(wrap_pyfunction!(calculate_distance_matrix, m)?)?;
    Ok(())
}