use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;
use std::time::Instant;

// dataset struct
struct DatasetData {
    gene_path: String,
    gem_path: String,
    gem_contains_cpg: bool
}

fn main() -> PyResult<()> {

    pyo3::prepare_freethreaded_python();

    // datasets
    let _small_dataset = DatasetData {
        gene_path: "tests/small_files/mRNA.csv".to_string(),
        gem_path: "tests/small_files/miRNA.csv".to_string(),
        gem_contains_cpg: false
    };

    let _methylation_dataset = DatasetData {
        gene_path: "tests/medium_files/methylation_gene.csv".to_string(),
        gem_path: "tests/medium_files/methylation_gem.csv".to_string(),
        gem_contains_cpg: true
    };

    let dataset_chosen = _methylation_dataset;

    let r1 = do_analysis(&dataset_chosen, AdjustmentMethod::BenjaminiHochberg);
    let r2 = do_analysis(&dataset_chosen, AdjustmentMethod::BenjaminiYekutieli);
    let r3 = do_analysis(&dataset_chosen, AdjustmentMethod::Bonferroni);

    if r1.is_ok() && r2.is_ok() && r3.is_ok() { Ok(()) } else { r1 } 
}

fn do_analysis(dataset_chosen: &DatasetData, adj_method: AdjustmentMethod) -> PyResult<()> {
    // Datasets's paths
    let gene_file_path = dataset_chosen.gene_path.clone();
    let gem_file_path = dataset_chosen.gem_path.clone();

    // Some parameters
    let gem_contains_cpg = dataset_chosen.gem_contains_cpg;       
    let is_all_vs_all = true;
    let keep_top_n = Some(10);
    let collect_gem_dataset = Some(true);  // false = disk (slow), true = ram, or none

    // Creates and run an analysis
    let analysis = Analysis {
        gene_file_path,
        gem_file_path,
        gem_contains_cpg,
        correlation_method: CorrelationMethod::Kendall,
        correlation_threshold: 0.6,
        sort_buf_size: 2_000_000,
        adjustment_method: adj_method,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n
    };

    let now = Instant::now();

    let (mut result, mut _total_combinations_count, mut number_of_combinations_evaluated) = (vec![], 0, 0);

    // repeat 10 times
    for _ in 0..10 {
        (result, _total_combinations_count, number_of_combinations_evaluated) =
        analysis.compute()?;
    }

    let milliseconds = now.elapsed().as_millis();

    /*
    if result.len() < 4 {
        for cor_p_value in result.iter() {
            println!("{}", cor_p_value);
        }
    } */

    println!("Finished in -> {} ms; average {} ms", milliseconds, milliseconds / 10);

    println!(
        "Number of elements -> {} of {} combinations evaluated",
        result.len(),
        number_of_combinations_evaluated
    );

    Ok(())
}