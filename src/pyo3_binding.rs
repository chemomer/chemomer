//! The Python binding (feature `pyo3`): the engine as a native
//! extension module `chemomer` with `from_molfile(molfile) ->
//! (inchi, inchikey)`.

// The `#[pyfunction]` wrapper's return conversion reads as an identity
// `PyErr` conversion to the current clippy (pyo3 0.22 macro output) —
// silenced here: macro-generated code outside this module's control.
#![allow(clippy::useless_conversion)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Derives `(inchi, inchikey)` from a V2000 molfile.
#[pyfunction]
fn from_molfile(molfile: &str) -> PyResult<(String, String)> {
    crate::Engine::from_molfile(molfile)
        .and_then(|engine| engine.inchi_string().map(|inchi| (engine, inchi)))
        .and_then(|(engine, inchi)| engine.inchikey().map(|inchikey| (inchi, inchikey)))
        .map_err(|error| PyValueError::new_err(error.message().to_string()))
}

/// Derives the AuxInfo line from a V2000 molfile.
#[pyfunction]
fn aux_from_molfile(molfile: &str) -> PyResult<String> {
    crate::Engine::from_molfile(molfile)
        .and_then(|engine| engine.aux_string())
        .map_err(|error| PyValueError::new_err(error.message().to_string()))
}

#[pymodule]
fn chemomer(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(from_molfile, module)?)?;
    module.add_function(wrap_pyfunction!(aux_from_molfile, module)?)?;
    Ok(())
}

// No in-crate unit tests here: a pyo3 test binary links the host
// libpython (unstable across interpreter installs), so the surface is
// verified where it runs — against the built extension module
// (scripts/smoke_py.py via maturin).
