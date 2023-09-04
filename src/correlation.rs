use bincode::{deserialize, serialize};
use extsort::Sortable;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyTuple};
use pyo3::ToPyObject;
use rgsl::randist::gaussian::gaussian_P;
use rgsl::{
    randist::t_distribution::{tdist_P, tdist_Q},
    statistics::{correlation, spearman},
};
use serde_derive::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::{
    fmt::Debug,
    io::{Read, Write},
};

/// Represents an correlation analysis result. Includes Gene, GEM, CpG Site ID (if specified) correlation statistic,
/// p-value and adjusted p-value.
#[pyclass(module = "ggca")]
#[derive(Clone, PartialEq, Serialize, Deserialize, Debug)]
pub struct CorResult {
    /// Gene name
    #[pyo3(get, set)]
    pub gene: String,
    /// Gene Expression Modulator (GEM) name
    #[pyo3(get, set)]
    pub gem: String,
    #[pyo3(get, set)]
    /// CpG Site ID
    pub cpg_site_id: Option<String>,
    /// Correlation statistic (Pearson, Spearman or Kendall, as selected)
    #[pyo3(get, set)]
    pub correlation: Option<f64>,
    /// P-value
    #[pyo3(get, set)]
    pub p_value: Option<f64>,
    /// Adjusted p-value (Benjamini-Hochberg, Benjamini-Yekutieli or Bonferroni, as selected)
    #[pyo3(get, set)]
    pub adjusted_p_value: Option<f64>,
}

#[pymethods]
impl CorResult {
    #[new]
    #[args(args = "*")]
    fn new(args: &PyTuple) -> Self {
        if args.len() >= 2 {
            CorResult {
                gene: args.get_item(0).unwrap().extract::<String>().unwrap(),
                gem: args.get_item(1).unwrap().extract::<String>().unwrap(),
                cpg_site_id: args
                    .get_item(2)
                    .unwrap()
                    .extract::<Option<String>>()
                    .unwrap(),
                correlation: args.get_item(3).unwrap().extract::<Option<f64>>().unwrap(),
                p_value: args.get_item(4).unwrap().extract::<Option<f64>>().unwrap(),
                adjusted_p_value: args.get_item(5).unwrap().extract::<Option<f64>>().unwrap(),
            }
        } else {
            CorResult {
                gene: String::from(""),
                gem: String::from(""),
                cpg_site_id: None,
                correlation: None,
                p_value: None,
                adjusted_p_value: None,
            }
        }
    }

    // Adds support for pickle
    pub fn __setstate__(&mut self, py: Python, state: PyObject) -> PyResult<()> {
        match state.extract::<&PyTuple>(py) {
            Ok(args) => {
                let gene_bytes = args.get_item(0).unwrap().extract::<&PyBytes>().unwrap();
                self.gene = deserialize(gene_bytes.as_bytes()).unwrap();
                let gem_bytes = args.get_item(1).unwrap().extract::<&PyBytes>().unwrap();
                self.gem = deserialize(gem_bytes.as_bytes()).unwrap();
                let cpg_site_id_bytes = args.get_item(2).unwrap().extract::<&PyBytes>().unwrap();
                self.cpg_site_id = deserialize(cpg_site_id_bytes.as_bytes()).unwrap();
                let correlation_bytes = args.get_item(3).unwrap().extract::<&PyBytes>().unwrap();
                self.correlation = deserialize(correlation_bytes.as_bytes()).unwrap();
                let p_value_bytes = args.get_item(4).unwrap().extract::<&PyBytes>().unwrap();
                self.p_value = deserialize(p_value_bytes.as_bytes()).unwrap();
                let adjusted_p_value_bytes =
                    args.get_item(5).unwrap().extract::<&PyBytes>().unwrap();
                self.adjusted_p_value = deserialize(adjusted_p_value_bytes.as_bytes()).unwrap();
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    // Adds support for pickle
    pub fn __getstate__(&self, py: Python) -> PyResult<PyObject> {
        let obj = (
            PyBytes::new(py, &serialize(&self.gene).unwrap()),
            PyBytes::new(py, &serialize(&self.gem).unwrap()),
            PyBytes::new(py, &serialize(&self.cpg_site_id).unwrap()),
            PyBytes::new(py, &serialize(&self.correlation).unwrap()),
            PyBytes::new(py, &serialize(&self.p_value).unwrap()),
            PyBytes::new(py, &serialize(&self.adjusted_p_value).unwrap()),
        )
            .to_object(py);
        Ok(obj)
    }

    /// Returns the CpG Site ID description. Empty String if it's None
    fn cpg_site_id_description(&self) -> &str {
        match &self.cpg_site_id {
            Some(cpg) => cpg,
            None => "",
        }
    }

    /// Gets the absolute correlation. Panics if the correlation value is None
    pub fn abs_correlation(&self) -> f64 {
        self.correlation.unwrap().abs()
    }

    // Will be auto-generated by PyO3
    pub fn __str__(&self) -> String {
        format!("{}", self)
    }

    // Will be auto-generated by PyO3
    pub fn __repr__(&self) -> String {
        format!(
            r#"CorResult("{}", "{}", "{}", {}, {:+e}, {:+e})"#,
            self.gene,
            self.gem,
            self.cpg_site_id_description(),
            self.correlation.unwrap_or(0.0),
            self.p_value.unwrap_or(0.0),
            self.adjusted_p_value.unwrap_or(0.0)
        )
    }
}

impl std::fmt::Display for CorResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            r#"Gene: "{}" | GEM: "{}" | CpG Site ID: "{}"
    Cor: {}
    P-value: {:+e}
    Adjusted p-value: {:+e}"#,
            self.gene,
            self.gem,
            self.cpg_site_id_description(),
            self.correlation.unwrap_or(0.0),
            self.p_value.unwrap_or(0.0),
            self.adjusted_p_value.unwrap_or(0.0)
        )
    }
}

impl Eq for CorResult {}

impl Sortable for CorResult {
    fn encode<W: Write>(&self, writer: &mut W) {
        let serialized = bincode::serialize(self).unwrap();
        writer.write_all(&serialized[..]).unwrap();
    }

    fn decode<R: Read>(reader: &mut R) -> Option<Self> {
        bincode::deserialize_from(reader).ok()
    }
}

pub trait Correlation: Sync {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64);
}

pub struct Pearson {
    n: usize,
    degrees_of_freedom: f64,
}

impl Pearson {
    fn new(n: usize) -> Self {
        Pearson {
            n,
            degrees_of_freedom: (n - 2) as f64,
        }
    }
}

impl Correlation for Pearson {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let r = correlation(x, 1, y, 1, self.n);

        // P-value (two-sided)
        // Based on R's cor.test method (https://github.com/SurajGupta/r-source/blob/a28e609e72ed7c47f6ddfbb86c85279a0750f0b7/src/library/stats/R/cor.test.R#L21)
        let statistic = self.degrees_of_freedom.sqrt() * r / (1.0 - r.powi(2)).sqrt();
        let p_value = 2.0
            * tdist_P(statistic, self.degrees_of_freedom)
                .min(tdist_Q(statistic, self.degrees_of_freedom));

        (r, p_value)
    }
}

struct Spearman {
    n: usize,
    degrees_of_freedom: f64,
}

impl Spearman {
    fn new(n: usize) -> Self {
        Spearman {
            n,
            degrees_of_freedom: (n - 2) as f64,
        }
    }
}

impl Correlation for Spearman {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let mut vec = Vec::with_capacity(2 * self.n);
        let workspace: &mut [f64] = vec.as_mut_slice();
        let rs = spearman(x, 1, y, 1, self.n, workspace);

        // P-value (two-sided)
        // Same behavior as Python Scipy's spearmanr method
        let t = rs * (self.degrees_of_freedom / ((rs + 1.0) * (1.0 - rs))).sqrt();
        let ccdf = tdist_Q(t.abs(), self.degrees_of_freedom);
        let p_value = 2.0 * ccdf;

        (rs, p_value)
    }
}

struct Kendall {}

impl Kendall {
    fn new(_n: usize) -> Self {
        Kendall {}
    }
}

impl Correlation for Kendall {
    fn correlate(&self, x: &[f64], y: &[f64]) -> (f64, f64) {
        let (tau, significance) = kendalls::tau_b_with_comparator(x, y, |a: &f64, b: &f64| {
            a.partial_cmp(b).unwrap_or(Ordering::Greater)
        })
        .unwrap();

        // P-value (two-sided)
        let cdf = gaussian_P(-significance.abs(), 1.0);
        let p_value = 2.0 * cdf;

        (tau, p_value)
    }
}

#[derive(Clone, Debug)]
pub enum CorrelationMethod {
    Spearman = 1,
    Kendall = 2,
    Pearson = 3,
}

impl std::fmt::Display for CorrelationMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let description = match &self {
            CorrelationMethod::Spearman => "Spearman",
            CorrelationMethod::Kendall => "Kendall",
            CorrelationMethod::Pearson => "Pearson",
        };

        write!(f, "{description}")
    }
}

pub fn get_correlation_method(
    correlation_method: &CorrelationMethod,
    number_of_samples: usize,
) -> Box<dyn Correlation> {
    match correlation_method {
        CorrelationMethod::Pearson => Box::new(Pearson::new(number_of_samples)),
        CorrelationMethod::Spearman => Box::new(Spearman::new(number_of_samples)),
        CorrelationMethod::Kendall => Box::new(Kendall::new(number_of_samples)),
    }
}
