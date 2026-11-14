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

// -- the reaction-identity surface (the RInChI layer) -------------------
// The RXNfile and RDfile input edges derive the oracle lines; the
// decode direction reads a RInChI string back into its components as
// the component table (`RinchiReaction::component_table`'s documented
// wire form: the flag line, TAB-separated `inchi<tab>aux<tab>key`
// triples, bare `<>` role separators).

/// Derives the RInChI identifier from a V2000 RXNfile.
#[wasm_bindgen]
pub fn rinchi_from_rxn(rxnfile: &str) -> Result<String, JsValue> {
    crate::parse_rxnfile(rxnfile)
        .and_then(|reaction| reaction.rinchi_string())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the RAuxInfo line from a V2000 RXNfile.
#[wasm_bindgen]
pub fn rinchi_aux_from_rxn(rxnfile: &str) -> Result<String, JsValue> {
    crate::parse_rxnfile(rxnfile)
        .and_then(|reaction| reaction.rinchi_auxinfo())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the short (hashed) RInChIKey from a V2000 RXNfile.
#[wasm_bindgen]
pub fn rinchi_short_key_from_rxn(rxnfile: &str) -> Result<String, JsValue> {
    crate::parse_rxnfile(rxnfile)
        .and_then(|reaction| reaction.rinchi_short_key())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the long RInChIKey from a V2000 RXNfile.
#[wasm_bindgen]
pub fn rinchi_long_key_from_rxn(rxnfile: &str) -> Result<String, JsValue> {
    crate::parse_rxnfile(rxnfile)
        .and_then(|reaction| reaction.rinchi_long_key())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the web RInChIKey from a V2000 RXNfile.
#[wasm_bindgen]
pub fn rinchi_web_key_from_rxn(rxnfile: &str) -> Result<String, JsValue> {
    crate::parse_rxnfile(rxnfile)
        .and_then(|reaction| reaction.rinchi_web_key())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Derives the RInChI identifier from an MDL RDfile (the `$DATUM
/// $MFMT` agent-scan input edge).
#[wasm_bindgen]
pub fn rinchi_from_rd(rdfile: &str) -> Result<String, JsValue> {
    crate::parse_rdfile(rdfile)
        .and_then(|reaction| reaction.rinchi_string())
        .map_err(|error| JsValue::from_str(error.message()))
}

/// Decodes a RInChI string (with its RAuxInfo companion; pass the
/// empty string when unknown) into the component table.
#[wasm_bindgen]
pub fn rinchi_decode(rinchi: &str, rauxinfo: &str) -> Result<String, JsValue> {
    crate::rinchi::reaction_from_rinchi_with_auxinfo(rinchi, rauxinfo)
        .and_then(|reaction| reaction.component_table())
        .map_err(|error| JsValue::from_str(error.message()))
}
