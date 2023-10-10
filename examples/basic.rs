use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;
use std::time::Instant;

fn main() -> PyResult<()> {
    pyo3::prepare_freethreaded_python();

    // Datasets's paths
    let gene_file_path = "tests/medium_files/methylation_gene.csv".to_string();
    let gem_file_path = "tests/medium_files/methylation_gem.csv".to_string();

    // Some parameters
    let gem_contains_cpg = true;
    let is_all_vs_all = true;
    let keep_top_n = Some(10); // Keeps the top 10 of correlation (sorting by abs values)
    let collect_gem_dataset = None; // Better performance. Keep small GEM files in memory

    let now = Instant::now();

    // Creates and run an analysis
    let analysis = Analysis {
        gene_file_path,
        gem_file_path,
        gem_contains_cpg,
        correlation_method: CorrelationMethod::Spearman,
        correlation_threshold: 0.7,
        sort_buf_size: 2_000_000,
        adjustment_method: AdjustmentMethod::BenjaminiHochberg,
        is_all_vs_all,
        collect_gem_dataset,
        keep_top_n,
    };

    let (result, _total_combinations_count, number_of_combinations_evaluated) =
        analysis.compute()?;

    let seconds = now.elapsed().as_secs();

    for cor_p_value in result.iter() {
        println!("{}", cor_p_value);
    }

    println!("Finished in -> {} seconds", seconds);
    println!(
        "Number of elements -> {} of {} combinations evaluated",
        result.len(),
        number_of_combinations_evaluated
    );

    Ok(())
}
