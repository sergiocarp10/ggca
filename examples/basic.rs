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
    let small_dataset = DatasetData {
        gene_path: "tests/small_files/mRNA.csv".to_string(),
        gem_path: "tests/small_files/miRNA.csv".to_string(),
        gem_contains_cpg: false
    };

    let _methylation_dataset = DatasetData {
        gene_path: "tests/medium_files/methylation_gene.csv".to_string(),
        gem_path: "tests/medium_files/methylation_gem.csv".to_string(),
        gem_contains_cpg: true
    };

    let dataset_chosen = small_dataset;

    // Datasets's paths
    let gene_file_path = dataset_chosen.gene_path;
    let gem_file_path = dataset_chosen.gem_path;

    // Some parameters
    let gem_contains_cpg = dataset_chosen.gem_contains_cpg;       
    let is_all_vs_all = true;
    let keep_top_n = None;
    let collect_gem_dataset = Some(true);  // false = disk (slow), true = ram, or none

    let now = Instant::now();

    // Creates and run an analysis
    let analysis = Analysis {
        gene_file_path,
        gem_file_path,
        gem_contains_cpg,
        correlation_method: CorrelationMethod::Pearson,
        correlation_threshold: 0.0,
        sort_buf_size: 2_000_000,
        adjustment_method: AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    };

    let (result, _total_combinations_count, number_of_combinations_evaluated) =
        analysis.compute()?;

    let milliseconds = now.elapsed().as_millis();

    // print results only if it isn't too much
    if result.len() < 4 {
        for cor_p_value in result.iter() {
            println!("{}", cor_p_value);
        }
    }

    println!("Finished in -> {} ms", milliseconds);

    println!(
        "Number of elements -> {} of {} combinations evaluated",
        result.len(),
        number_of_combinations_evaluated
    );

    Ok(())
}
