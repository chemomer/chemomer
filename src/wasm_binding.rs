//! The WebAssembly binding (feature `wasm`): the engine for the site
//! playground and browser consumers.

use wasm_bindgen::prelude::*;

/// Derives the identifier from a V2000 molfile.
#[wasm_bindgen]
pub fn from_molfile(molfile: &str) -> Result<String, JsValue> {
    crate::Engine::from_molfile(molfile)
        .and_then(|engine| engine.inchi_string())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the AuxInfo line from a V2000 molfile.
#[wasm_bindgen]
pub fn aux_from_molfile(molfile: &str) -> Result<String, JsValue> {
    crate::Engine::from_molfile(molfile)
        .and_then(|engine| engine.aux_string())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the InChIKey from a V2000 molfile.
#[wasm_bindgen]
pub fn inchikey(molfile: &str) -> Result<String, JsValue> {
    crate::Engine::from_molfile(molfile)
        .and_then(|engine| engine.inchikey())
        .map_err(|error| JsValue::from_str(error.message()))
}
