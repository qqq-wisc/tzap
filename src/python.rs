//! PyO3 entry points for the mixed Rust/Python package.
//!
//! The public, typed API lives in `python/tzap`. This module deliberately
//! exposes a small private boundary: QASM in, QASM plus primitive metrics out.
//! Keeping Python-specific types here means the optimizer and its Rust API
//! remain completely independent of Python.

use pyo3::create_exception;
use pyo3::exceptions::{PyException, PyValueError};
use pyo3::prelude::*;

use crate::circuit::Circuit;
use crate::optimize::{
    Level, Metrics, Options, PassName, Report, SuperOptBounds, SuperOptGates, optimize,
};

create_exception!(_native, TzapError, PyException);
create_exception!(_native, QasmError, TzapError);
create_exception!(_native, OptimizationError, TzapError);

type RawMetrics = (usize, usize, usize, usize, usize);
type RawReport = (RawMetrics, RawMetrics, RawMetrics);

fn raw_metrics(metrics: Metrics) -> RawMetrics {
    (
        metrics.gates,
        metrics.two_qubit,
        metrics.depth,
        metrics.t,
        metrics.rz,
    )
}

fn raw_report(report: Report) -> RawReport {
    (
        raw_metrics(report.input),
        raw_metrics(report.baseline),
        raw_metrics(report.output),
    )
}

fn parse_level(level: &str) -> PyResult<Level> {
    match level.to_ascii_lowercase().as_str() {
        "o1" | "1" => Ok(Level::O1),
        "o2" | "2" => Ok(Level::O2),
        "o3" | "3" => Ok(Level::O3),
        "osuper" | "super" => Ok(Level::Osuper),
        _ => Err(PyValueError::new_err(format!(
            "unknown optimization level {level:?}; expected 'O1', 'O2', 'O3', or 'Osuper'"
        ))),
    }
}

fn parse_passes(passes: Option<Vec<String>>) -> PyResult<Option<Vec<PassName>>> {
    passes
        .map(|passes| {
            if passes.is_empty() {
                return Err(PyValueError::new_err(
                    "passes must contain at least one pass name",
                ));
            }
            passes
                .into_iter()
                .map(|name| {
                    PassName::parse(&name).ok_or_else(|| {
                        PyValueError::new_err(format!(
                            "unknown pass {name:?}; available passes: {}",
                            PassName::all_names()
                        ))
                    })
                })
                .collect()
        })
        .transpose()
}

fn positive_bound(value: Option<usize>, name: &str) -> PyResult<Option<usize>> {
    if value == Some(0) {
        Err(PyValueError::new_err(format!(
            "{name} must be a positive integer"
        )))
    } else {
        Ok(value)
    }
}

/// Run the same optimizer driver as the CLI.
///
/// This is private because `python/tzap/_core.py` wraps the primitive report
/// tuple in stable Python dataclasses and documents the public contract.
#[pyfunction]
#[pyo3(signature = (
    qasm,
    *,
    level = "O3",
    passes = None,
    fixpoint = false,
    decompose_rz = false,
    decompose_cz = false,
    decompose_ccx = false,
    rz_epsilon = crate::optimize::DEFAULT_RZ_EPSILON,
    parallel = false,
    superopt_qubits = None,
    superopt_window_gates = None,
    superopt_murm_entries = None,
    superopt_gates = "auto",
))]
#[allow(clippy::too_many_arguments)]
fn _optimize_qasm(
    py: Python<'_>,
    qasm: &str,
    level: &str,
    passes: Option<Vec<String>>,
    fixpoint: bool,
    decompose_rz: bool,
    decompose_cz: bool,
    decompose_ccx: bool,
    rz_epsilon: f64,
    parallel: bool,
    superopt_qubits: Option<usize>,
    superopt_window_gates: Option<usize>,
    superopt_murm_entries: Option<usize>,
    superopt_gates: &str,
) -> PyResult<(String, RawReport)> {
    if !rz_epsilon.is_finite() || rz_epsilon <= 0.0 {
        return Err(PyValueError::new_err(
            "rz_epsilon must be a positive, finite number",
        ));
    }

    let passes = parse_passes(passes)?;
    if passes.is_some() && (decompose_rz || decompose_cz || decompose_ccx) {
        return Err(PyValueError::new_err(
            "passes cannot be combined with decomposition options; include the decomposition passes instead",
        ));
    }

    let options = Options {
        level: parse_level(level)?,
        passes,
        fixpoint,
        decompose_rz,
        decompose_cz,
        decompose_ccx,
        rz_epsilon,
        parallel,
        superopt: SuperOptBounds {
            qubits: positive_bound(superopt_qubits, "superopt_qubits")?,
            window_gates: positive_bound(superopt_window_gates, "superopt_window_gates")?,
            murm_entries: positive_bound(superopt_murm_entries, "superopt_murm_entries")?,
        },
        superopt_gates: SuperOptGates::parse(superopt_gates).map_err(PyValueError::new_err)?,
    };

    let circuit = Circuit::from_qasm(qasm).map_err(QasmError::new_err)?;
    let (optimized, report) = py
        .detach(|| optimize(&circuit, &options))
        .map_err(|error| OptimizationError::new_err(error.to_string()))?;
    Ok((optimized.to_qasm(), raw_report(report)))
}

#[pymodule]
#[pyo3(name = "_native")]
fn python_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_optimize_qasm, module)?)?;
    let tzap_error = module.py().get_type::<TzapError>();
    let qasm_error = module.py().get_type::<QasmError>();
    let optimization_error = module.py().get_type::<OptimizationError>();
    for exception in [&tzap_error, &qasm_error, &optimization_error] {
        exception.setattr("__module__", "tzap._native")?;
    }
    module.add("TzapError", tzap_error)?;
    module.add("QasmError", qasm_error)?;
    module.add("OptimizationError", optimization_error)?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
