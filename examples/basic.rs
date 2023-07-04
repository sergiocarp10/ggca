use ggca::adjustment::AdjustmentMethod;
use ggca::analysis::Analysis;
use ggca::correlation::CorrelationMethod;
use pyo3::PyResult;
use std::time::Instant;

// Datasets's paths
const GENE_SMALL_FILE_PATH: &str = "tests/small_files/mRNA.csv"; // mRNA = 600 rows
const GEM_SMALL_FILE_PATH: &str = "tests/small_files/miRNA.csv"; // miRNA = 299 rows

const GENE_METHYLATION_FILE_PATH: &str = "tests/medium_files/methylation_gene.csv"; // mRNA = 41 rows
const METHYLATION_FILE_PATH: &str = "tests/medium_files/methylation_gem.csv"; // miRNA = 1505 rows


fn main() -> PyResult<()> {

    pyo3::prepare_freethreaded_python();

    // Datasets's paths
    let gene_file_path = GENE_SMALL_FILE_PATH.to_string();
    let gem_file_path = GEM_SMALL_FILE_PATH.to_string();

    // Some parameters
    let gem_contains_cpg = false;       // false en small, true en methylation
    let is_all_vs_all = true;
    let keep_top_n = None;
    let collect_gem_dataset = Some(false);  // false = disk (slow), true = ram, or none

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
