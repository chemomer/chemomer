//! chemomer — the chemical molecule model + identity engine,
//! gated byte-for-byte by the oracle corpus in
//! [inchi-tests](https://github.com/chemomer/inchi-tests).
//!
//! Architecture: an [`Engine`] façade over encapsulated concerns —
//! `molfile` reads the input edge, `graph` owns molecular-graph
//! semantics, `canonical` owns canonical ranking, and `layers` holds
//! one emitter per InChI layer (open for extension, closed for
//! modification).

// Unsafe code is denied everywhere except the C ABI (src/ffi.rs),
// where it is the job — the boundary module carries its own allow.
#![deny(unsafe_code)]

#[allow(unsafe_code)]
pub mod ffi;

#[cfg(feature = "napi")]
pub mod napi_binding;
#[cfg(feature = "pyo3")]
pub mod pyo3_binding;
#[cfg(feature = "wasm")]
pub mod wasm_binding;

pub mod aux;
pub mod canonical;
pub mod engine;
mod error;
pub mod graph;
pub mod key;
pub mod layers;
mod molfile;
pub mod options;
pub mod rinchi;

pub use engine::Engine;
pub use error::{Error, Result};
pub use graph::{AtomRecord, Structure};
pub use options::Options;
pub use rinchi::{parse_rxn, RinchiComponent, RinchiReaction};

/// Parses the V2000 molfile subset into a [`Structure`].
pub fn parse_molfile(input: &str) -> Result<Structure> {
    molfile::parse(input)
}

/// Parses an MDL V2000 RXNfile into a [`RinchiReaction`] — the
/// reaction-identity layer's input edge.
pub fn parse_rxnfile(input: &str) -> Result<RinchiReaction> {
    rinchi::parse_rxn(input)
}
